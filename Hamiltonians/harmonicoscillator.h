#ifndef __HARMONIC_OSCILLATOR__
#define __HARMONIC_OSCILLATOR__
#include <memory>
#include <vector>

#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(double omega);
    double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
    );

private:
    double m_omega;
};
#endif