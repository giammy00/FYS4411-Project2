#include <iostream>
#include <vector>
#include <memory>
#include <math.h>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

using namespace std;

/*
h-bar = 1
m = 1
omega = 1
*/
std::unique_ptr<Sampler> runSimulation(
    unsigned int numberOfDimensions,
    unsigned int numberOfParticles,
    unsigned int numberOfMetropolisSteps,
    unsigned int numberOfEquilibrationSteps,
    double omega,
    double a_ho,
    double alpha,
    double stepLength
){
    int seed = 2023;
    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng, a_ho);//[x]
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(//[x]
            // Construct unique_ptr to Hamiltonian
            std::make_unique<HarmonicOscillator>(omega),//[x]
            // Construct unique_ptr to wave function
            std::make_unique<SimpleGaussian>(alpha),//[x]
            // Construct unique_ptr to solver, and move rng
            std::make_unique<Metropolis>(std::move(rng)),//[x]
            // Move the vector of particles to system
            std::move(particles));

    // Run steps to equilibrate particles
    auto sampler = system->runEquilibrationSteps(//[x]
            stepLength,
            numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    sampler = system->runMetropolisSteps(//[ ]
            std::move(sampler),
            stepLength,
            numberOfMetropolisSteps);

    // Return the results
    return sampler;
}

int main() {
    // Seed for the random number generator
    int seed = 2023;

    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 1;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e6;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e5;
    double omega = 1.0; // Oscillator frequency.
    double a_ho = std::sqrt(1./omega); // Characteristic size of the Harmonic Oscillator
    double alpha = 0.5; // Variational parameter.
    double stepLength = 0.5; // Metropolis step length.
    stepLength *= a_ho; // Scale the steplength in case of changed omega

    for(alpha = alpha; alpha < 2; alpha += 0.05){

        auto sampler = runSimulation(
                numberOfDimensions,
                numberOfParticles,
                numberOfMetropolisSteps,
                numberOfEquilibrationSteps,
                omega,
                a_ho,
                alpha,
                stepLength);

        // Output information from the simulation
        sampler->printOutputToTerminalShort();//[ ]
    }
    return 0;
}
