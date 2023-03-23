#include <memory>
#include <vector>

#include "metropolis.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"

#include <iostream>
using std::cout;
using std::endl;

Metropolis::Metropolis(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool Metropolis::step(
        double stepLength,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles)
{
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change its position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated at
     * this new position with the one at the old position).
     */
    int index = m_rng->nextInt(0,particles.size()-1);
    
    auto step = std::vector<double>();
    for(unsigned int i=0; i<particles[0]->getNumberOfDimensions(); i++){
        step.push_back(m_rng->nextGaussian(0,stepLength));
    }
    double q = waveFunction.phiRatio(particles, index, step);


    if(m_rng->nextDouble()<q){
        particles[index]->adjustPosition(step);
        waveFunction.updateCachedVariables(step);
        return true;
    }
    
    return false;
}
