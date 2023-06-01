#include<memory>
#include <cassert>
#include <iostream>

#include "harmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

// No interaction
// Spherical geometry

HarmonicOscillator::HarmonicOscillator(double omega)
{
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        )
{
    /* Here, you need to compute the kinetic and potential energies.
     * Access to the wave function methods can be done using the dot notation
     * for references, e.g., wavefunction.computeDoubleDerivative(particles),
     * to get the Laplacian of the wave function.
     * */

    double r2 = 0;
    unsigned int N = particles.size();
    for (unsigned int i = 0; i < N; i++){
        auto position = particles[i]->getPosition();
        for (unsigned int i = 0; i<particles[0]->getNumberOfDimensions(); i++)
            r2 += position[i]*position[i];
    }
    // m = omega = 1
    double potentialEnergy = 0.5 * r2;
    double kineticEnergy   = waveFunction.computeDoubleDerivative(particles)*-0.5;
    return (kineticEnergy + potentialEnergy)/N;
}
