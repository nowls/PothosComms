// Copyright (c) 2014-2016 Josh Blum
// SPDX-License-Identifier: BSL-1.0

#include <Pothos/Framework.hpp>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <complex>
#include <random>

static const size_t waveTableSize = 4096;

/***********************************************************************
 * |PothosDoc Noise Source
 *
 * The noise source produces pseudorandom noise patterns.
 * When a complex data type is chosen, the real and imaginary
 * components are simply treated as two independent channels.
 *
 * |category /Sources
 * |category /Waveforms
 * |category /Random
 * |keywords noise random source pseudorandom gaussian
 * |alias /blocks/noise_source
 *
 * |param dtype[Data Type] The data type produced by the noise source.
 * |widget DTypeChooser(float=1,cfloat=1,int=1,cint=1)
 * |default "complex_float32"
 * |preview disable
 *
 * |param wave[Wave Type] The type of the pseudorandom noise produced.
 * |option [Uniform] "UNIFORM"
 * |option [Normal] "NORMAL"
 * |option [Laplace] "LAPLACE"
 * |option [Poisson] "POISSON"
 * |default "NORMAL"
 *
 * |param ampl[Amplitude] A constant scalar representing the amplitude.
 * |default 1.0
 *
 * |param offset A constant value added to the waveform after scaling.
 * |default 0.0
 * |preview valid
 *
 * |param mean The mean of the distribution - applies to all distributions.
 * |default 0.0
 * |preview valid
 *
 * |param b A value with distribution-dependent meaning:
 * <ul>
 *   <li><b>Uniform distribution:</b> range = mean +/- b</li>
 *   <li><b>Normal distribution:</b> the standard deviation</li>
 *   <li><b>Laplace distribution:</b> the diversity parameter</li>
 *   <li><b>Poisson distribution:</b> not used</li>
 * </ul>
 * |default 1.0
 *
 * |param fast[Fast] Use a sample pool for speed at the expense of less entropy.
 * |option [Enabled] true
 * |option [Disabled] false
 * |default true
 * |preview invalid
 *
 * |factory /comms/noise_source(dtype)
 * |setter setWaveform(wave)
 * |setter setOffset(offset)
 * |setter setAmplitude(ampl)
 * |setter setMean(mean)
 * |setter setB(b)
 **********************************************************************/
template <typename Type>
class NoiseSource : public Pothos::Block
{
public:
    NoiseSource(void):
        _index(0),
        _table(waveTableSize),
        _offset(0.0),
        _scalar(1.0),
        _wave("NORMAL"),
        _mean(0.0),
        _b(1.0),
        _fast(true),
        _gen(_rd()),
        _waveIndex(0, waveTableSize-1)
    {
        this->setupOutput(0, typeid(Type));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, setWaveform));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, getWaveform));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, setOffset));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, getOffset));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, setAmplitude));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, getAmplitude));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, setMean));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, getMean));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, setB));
        this->registerCall(this, POTHOS_FCN_TUPLE(NoiseSource, getB));
    }

    void activate(void)
    {
        this->updateTable();
    }

    void work(void)
    {
        auto outPort = this->output(0);
        Type *out = outPort->buffer();
        if (_fast)
        {
            _index += _waveIndex(_gen); //lookup into table is random each work()
            for (size_t i = 0; i < outPort->elements(); i++)
            {
                out[i] = _table[_index % waveTableSize];
                _index++;
            }
        }
        else
        {
            for (size_t i = 0; i < outPort->elements(); i++)
            {
                if (_wave == "UNIFORM") this->setElem(out[i], std::complex<double>(_uniform(_gen), _uniform(_gen)));
                else if (_wave == "NORMAL") this->setElem(out[i], std::complex<double>(_normal(_gen), _normal(_gen)));
                else if (_wave == "LAPLACE") this->setElem(out[i], std::complex<double>(_laplace(_gen), _laplace(_gen)));
                else if (_wave == "POISSON") this->setElem(out[i], std::complex<double>(_poisson(_gen), _poisson(_gen)));
                else throw Pothos::Exception("NoiseSource::work(void)", "unknown waveform setting: "+_wave);
            }
        }
        outPort->produce(outPort->elements());
    }

    void setWaveform(const std::string &wave)
    {
        _wave = wave;
        this->updateTable();
    }

    std::string getWaveform(void)
    {
        return _wave;
    }

    void setOffset(const std::complex<double> &offset)
    {
        _offset = offset;
        this->updateTable();
    }

    std::complex<double> getOffset(void) const
    {
        return _offset;
    }

    void setAmplitude(const std::complex<double> &scalar)
    {
        _scalar = scalar;
        this->updateTable();
    }

    std::complex<double> getAmplitude(void) const
    {
        return _scalar;
    }

    void setMean(const double mean)
    {
        _mean = mean;
        this->updateTable();
    }

    double getMean(void) const
    {
        return _mean;
    }

    void setB(const double b)
    {
        _b = b;
        this->updateTable();
    }

    double getB(void) const
    {
        return _b;
    }

private:
    void updateTable(void)
    {
        if (not this->isActive()) return;

        if (_wave == "UNIFORM")
        {
            _uniform = std::uniform_real_distribution<>(_mean-_b, _mean+_b);
            for (size_t i = 0; i < _table.size(); i++)
            {
                this->setElem(_table[i], std::complex<double>(_uniform(_gen), _uniform(_gen)));
            }
        }
        else if (_wave == "NORMAL")
        {
            _normal = std::normal_distribution<>(_mean, _b);
            for (size_t i = 0; i < _table.size(); i++)
            {
                this->setElem(_table[i], std::complex<double>(_normal(_gen), _normal(_gen)));
            }
        }
        else if (_wave == "LAPLACE")
        {
            _uniform = std::uniform_real_distribution<>(_mean-_b, _mean+_b);
            for (size_t i = 0; i < _table.size(); i++)
            {
                this->setElem(_table[i], std::complex<double>(_laplace(_gen), _laplace(_gen)));
            }
        }
        else if (_wave == "POISSON")
        {
            _poisson = std::poisson_distribution<>(_mean);
            for (size_t i = 0; i < _table.size(); i++)
            {
                this->setElem(_table[i], std::complex<double>(_poisson(_gen), _poisson(_gen)));
            }
        }
        else throw Pothos::InvalidArgumentException("NoiseSource::setWaveform("+_wave+")", "unknown waveform setting");
    }

    template <typename T>
    void setElem(T &out, const std::complex<double> &val)
    {
        out = Type((_scalar * val + _offset).real());
    }

    template <typename T>
    void setElem(std::complex<T> &out, const std::complex<double> &val)
    {
        out = Type(_scalar * val + _offset);
    }

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template <typename GenType>
    double _laplace(GenType &gen)
    {
        //http://en.wikipedia.org/wiki/Laplace_distribution
        auto num = _uniform(gen);
        if (num < 0) return _mean + _b*std::log(1+num);
        else return _mean - _b*std::log(1-num);
    }

    size_t _index;
    std::vector<Type> _table;
    std::complex<double> _offset, _scalar;
    std::string _wave;
    double _mean;
    double _b;
    bool _fast;

    std::random_device _rd;
    std::mt19937 _gen;
    std::uniform_int_distribution<size_t> _waveIndex;
    std::uniform_real_distribution<> _uniform;
    std::normal_distribution<> _normal;
    std::poisson_distribution<> _poisson;
};

/***********************************************************************
 * registration
 **********************************************************************/
static Pothos::Block *noiseSourceFactory(const Pothos::DType &dtype)
{
    #define ifTypeDeclareFactory(type) \
        if (dtype == Pothos::DType(typeid(type))) return new NoiseSource<type>(); \
        if (dtype == Pothos::DType(typeid(std::complex<type>))) return new NoiseSource<std::complex<type>>();
    ifTypeDeclareFactory(double);
    ifTypeDeclareFactory(float);
    ifTypeDeclareFactory(int64_t);
    ifTypeDeclareFactory(int32_t);
    ifTypeDeclareFactory(int16_t);
    ifTypeDeclareFactory(int8_t);
    throw Pothos::InvalidArgumentException("noiseSourceFactory("+dtype.toString()+")", "unsupported type");
}

static Pothos::BlockRegistry registerNoiseSource(
    "/comms/noise_source", &noiseSourceFactory);

static Pothos::BlockRegistry registerNoiseSourceOldPath(
    "/blocks/noise_source", &noiseSourceFactory);
