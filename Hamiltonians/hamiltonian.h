#ifndef __HAMILTONIAN__
#define __HAMILTONIAN__
#include <memory>
#include <vector>

class Hamiltonian {
public:
    virtual ~Hamiltonian() = default;
    virtual double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
    ) = 0;
};

#endif