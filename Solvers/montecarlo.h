#ifndef __MONTECARLO__
#define __MONTECARLO__
#include <vector>
#include <memory>

class MonteCarlo {
public:
    MonteCarlo(std::unique_ptr<class Random> rng);
    virtual ~MonteCarlo() = default;

    virtual bool step( // Only placeholder, will always fail
            double stepLength,
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles) = 0;

protected:
    std::unique_ptr<class Random> m_rng;
};
#endif