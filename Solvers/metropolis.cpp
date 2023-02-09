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
    int dimension = m_rng->nextInt(0, particles[0]->getNumberOfDimensions()-1);
    double step = m_rng->nextGaussian(0., stepLength);
    // TODO make function to compare energies of adjacent states [ ]
    double old_phi = waveFunction.evaluate(particles);
    old_phi *= old_phi;
    particles[index]->adjustPosition(step, dimension);
    double new_phi = waveFunction.evaluate(particles);
    new_phi *= new_phi;

    // cout << "old_phi = " << old_phi << endl;
    // cout << "new_phi = " << new_phi << endl;
    // cout << "ratio = " << new_phi/old_phi << endl;

    if(new_phi>old_phi || m_rng->nextDouble()<(new_phi/old_phi)){
        // cout << "accepted" << endl;
        return true;
    }
    
    particles[index]->adjustPosition(-step, dimension);
    // cout << "denied" << endl;
    return false;
}
