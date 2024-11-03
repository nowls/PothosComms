// Copyright (c) 2014-2019 Josh Blum
// SPDX-License-Identifier: BSL-1.0

#include <Pothos/Framework.hpp>
#include <Pothos/Util/QFormat.hpp>
#include <cstdint>
#include <complex>
#include <iostream>
#include <algorithm> //min/max
#include <chrono>
#include <limits>

/***********************************************************************
 * |PothosDoc Signal Probe
 *
 * The signal probe block records the last calculation from a stream of elements.
 * The signal probe has a slot called "probeValue" will will cause
 * a signal named "valueTriggered" to emit the most recent value.
 * The probe will also emit the value automatically at the specified rate
 * using the "valueChanged" signal.
 *
 * The calculation for value can be, the last seen value,
 * the RMS (root mean square) over the last buffer,
 * the mean (average value) over the last buffer,
 * the minimum magnitude over the last buffer,
 * the maximum magnitude over the last buffer,
 * the P2RMS (peak-to-RMS ratio) over the last buffer,
 * or the SINAD (signal to noise-and-distortion ratio) over the last buffer.
 *
 * |category /Utility
 * |category /Event
 * |keywords rms average mean min minimum max maximum papr snr sinad
 * |alias /blocks/stream_probe
 *
 * |param dtype[Data Type] The data type consumed by the stream probe.
 * |widget DTypeChooser(float=1,cfloat=1,int=1,cint=1)
 * |default "complex_float32"
 * |preview disable
 *
 * |param mode The calculation mode for the value.
 * In value mode, this block expects to be fed by an upstream block
 * that produces a stream of slow-changing values.
 * Otherwise the value will appear random.
 * Note that the SINAD mode is most useful when the input contains samples from
 * an FFT block with its size equal to the signal probe window length.
 * |default "VALUE"
 * |option [Value] "VALUE"
 * |option [RMS] "RMS"
 * |option [Mean] "MEAN"
 * |option [Min] "MIN"
 * |option [Max] "MAX"
 * |option [P2RMS] "P2RMS"
 * |option [SINAD] "SINAD"
 *
 * |param rate How many calculations per second?
 * The probe will perform a calculation at most this many times per second.
 * Incoming samples will be dropped and not processed between calculations.
 * A special value of 0.0 means perform the calculation on every input window.
 * |preview valid
 * |default 0.0
 *
 * |param window How many elements to calculate over?
 * |default 1024
 *
 * |factory /comms/signal_probe(dtype)
 * |setter setMode(mode)
 * |setter setRate(rate)
 * |setter setWindow(window)
 **********************************************************************/
template <typename Type, typename ProbeType>
class SignalProbe : public Pothos::Block
{
public:
    SignalProbe(void):
        _value(0),
        _mode("VALUE"),
        _window(1024),
        _rate(0.0)
    {
        this->setupInput(0, typeid(Type));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, value));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, setMode));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, getMode));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, setWindow));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, getWindow));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, setRate));
        this->registerCall(this, POTHOS_FCN_TUPLE(SignalProbe, getRate));
        this->registerProbe("value");
        this->registerSignal("valueChanged");
        this->input(0)->setReserve(1);
    }

    ProbeType value(void)
    {
        return _value;
    }

    void setMode(const std::string &mode)
    {
        _mode = mode;
    }

    std::string getMode(void) const
    {
        return _mode;
    }

    void setWindow(const size_t window)
    {
        _window = window;
        this->input(0)->setReserve(window);
    }

    size_t getWindow(void) const
    {
        return _window;
    }

    void setRate(const double rate)
    {
        _rate = rate;
    }

    double getRate(void) const
    {
        return _rate;
    }

    void activate(void)
    {
        _nextCalc = std::chrono::high_resolution_clock::now();
    }

    void work(void)
    {
        auto inPort = this->input(0);
        const Type *x = inPort->buffer();
        const auto N = std::min(_window, inPort->elements());
        inPort->consume(N);

        //check if the time expired or rate is 0.0
        auto currentTime = std::chrono::high_resolution_clock::now();
        if (_rate != 0.0 and currentTime < _nextCalc) return;

        //increment for the next calculation time
        if (_rate != 0.0)
        {
            const auto tps = std::chrono::nanoseconds((long long)(1e9/_rate));
            _nextCalc += std::chrono::duration_cast<std::chrono::high_resolution_clock::duration>(tps);
        }

        if (_mode == "VALUE") _value = Pothos::Util::fromQ<ProbeType>(x[N-1], 0);
        else if (_mode == "RMS") _value = computeRMS(x, N);
        else if (_mode == "MEAN")
        {
            ProbeType mean = 0;
            for (size_t n = 0; n < N; n++) mean += Pothos::Util::fromQ<ProbeType>(x[n], 0);
            mean /= N;
            _value = mean;
        }
        else if (_mode == "MIN")
        {
            double minimum = std::numeric_limits<double>::max();
            ProbeType x_n;
            for (size_t n = 0; n < N; n++)
            {
                x_n = Pothos::Util::fromQ<ProbeType>(x[n], 0);
                const double v = std::abs(x_n);
                if (v < minimum) minimum = v;
            }
            _value = minimum;
        }
        else if (_mode == "MAX") _value = computeMaxMag(x, N);
        else if (_mode == "P2RMS") _value = computeMaxMag(x, N) / computeRMS(x, N);
        else if (_mode == "SINAD")
        {
            size_t index = computeMaxMagIndex(x, N);
            if (index == 0) index++;
            else if (index == N - 1) index--;
            double rms = std::abs(computeRMS(x, N));
            double rms_sig = std::abs(computeRMS(&x[index - 1], 3));

            double rss_nd = std::sqrt(double(N)*rms*rms - 3.0*rms_sig*rms_sig);
            _value = rms_sig / rss_nd;
        }

        this->emitSignal("valueChanged", _value);
    }

private:
    ProbeType computeRMS(const Type *x, const size_t N)
    {
        double accumulator = 0.0;
        ProbeType x_n;
        for (size_t n = 0; n < N; n++)
        {
            x_n = Pothos::Util::fromQ<ProbeType>(x[n], 0);
            const double v = std::abs(x_n);
            accumulator += v*v;
        }
        return std::sqrt(accumulator/N);
    }

    size_t computeMaxMagIndex(const Type *x, const size_t N)
    {
        double maximum = std::numeric_limits<double>::lowest();
        size_t index = 0;
        ProbeType x_n;
        for (size_t n = 0; n < N; n++)
        {
            x_n = Pothos::Util::fromQ<ProbeType>(x[n], 0);
            const double v = std::abs(x_n);
            if (v > maximum)
            {
                index = n;
                maximum = v;
            }
        }
        return index;
    }

    ProbeType computeMaxMag(const Type *x, const size_t N)
    {
        const size_t index = computeMaxMagIndex(x, N);
        return std::abs(Pothos::Util::fromQ<ProbeType>(x[index], 0));
    }

    ProbeType _value;
    std::string _mode;
    size_t _window;
    double _rate;
    std::chrono::high_resolution_clock::time_point _nextCalc;
};

/***********************************************************************
 * registration
 **********************************************************************/
static Pothos::Block *signalProbeFactory(const Pothos::DType &dtype)
{
    #define ifTypeDeclareFactory(type) \
        if (dtype == Pothos::DType(typeid(type))) return new SignalProbe<type, double>(); \
        if (dtype == Pothos::DType(typeid(std::complex<type>))) return new SignalProbe<std::complex<type>, std::complex<double>>();
    ifTypeDeclareFactory(double);
    ifTypeDeclareFactory(float);
    ifTypeDeclareFactory(int64_t);
    ifTypeDeclareFactory(int32_t);
    ifTypeDeclareFactory(int16_t);
    ifTypeDeclareFactory(int8_t);
    throw Pothos::InvalidArgumentException("signalProbeFactory("+dtype.toString()+")", "unsupported type");
}

static Pothos::BlockRegistry registerSignalProbe(
    "/comms/signal_probe", &signalProbeFactory);

static Pothos::BlockRegistry registerSignalProbeOldPath(
    "/blocks/stream_probe", &signalProbeFactory);
