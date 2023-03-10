#include <iostream>
#include <memory>
#include <cassert>

#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"


System::System(
        std::unique_ptr<class Hamiltonian> hamiltonian,
        std::unique_ptr<class WaveFunction> waveFunction,
        std::unique_ptr<class MonteCarlo> solver,
        std::vector<std::unique_ptr<class Particle>> particles)
{
    m_numberOfParticles = particles.size();;
    m_numberOfDimensions = particles[0]->getNumberOfDimensions();
    m_hamiltonian = std::move(hamiltonian);
    m_waveFunction = std::move(waveFunction);
    m_solver = std::move(solver);
    m_particles = std::move(particles);
}


std::unique_ptr<class Sampler> System::runEquilibrationSteps(
        double stepLength,
        unsigned int numberOfEquilibrationSteps)
{
    auto sampler = std::make_unique<Sampler>(
            m_numberOfParticles,
            m_numberOfDimensions, 
            m_waveFunction->getNumberOfParameters()
            );

    for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++) {
        sampler->equilibrationSample(m_solver->step(stepLength, *m_waveFunction, m_particles));
    }

    return sampler;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
        std::unique_ptr<class Sampler> sampler,
        double stepLength,
        unsigned int numberOfMetropolisSteps)
{
    for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
        /* Call solver method to do a single Monte-Carlo step.
         */
        bool acceptedStep = m_solver->step(stepLength, *m_waveFunction, m_particles);

        /* Here you should sample the energy (and maybe other things) using the
         * sampler instance of the Sampler class.
         */
        sampler->sample(acceptedStep, this);
    }

    sampler->transferWaveFunctionParameters(m_waveFunction->getParameters());
    sampler->computeAverages();

    return sampler;
}

double System::computeLocalEnergy()
{
    // Helper function
    return m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles);
}

const std::vector<double>& System::getWaveFunctionParameters()
{
    // Helper function
    return m_waveFunction->getParameters();
}


//Return sampled quantities for computation of gradient
std::vector<double> System::getdPhi_dParams()
{
    return m_waveFunction->getdPhi_dParams();
}

//compute gradient of local energy wrt a variational parameter
//nb: avgGradientTerms must contain  MC estimates of <dLogPsi/d(parameter)> AND <dLogPsi/d(parameter)*LocalEnergy>
double System::computeDerivative( double avgElocal, std::vector<double> avgGradientTerms)
{
    return 2*(avgGradientTerms[1]-avgElocal*avgGradientTerms[0]);
}