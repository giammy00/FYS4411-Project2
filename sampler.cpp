#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions
        )
{
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
}

void Sampler::equilibrationSample(bool acceptedStep){
    m_equilibrationStepNumber++;
    m_numberOfAcceptedEquilibrationSteps += acceptedStep;
}

void Sampler::sample(bool acceptedStep, System* system) {
    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += localEnergy * localEnergy;
    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;

    //HERE WE NEED TO CARE ABOUT **WHICH** QUANTITIES TO SAMPLE 
    //In fact, they  dependon the wavefunction used.
    //system->waveFunction.computeGradientTerms()
}

void Sampler::transferWaveFunctionParameters(std::vector<double> parameters){
    m_waveFunctionParameters = parameters;
}

void Sampler::printOutputToTerminal() {
    auto pa = m_waveFunctionParameters;
    auto p = pa.size();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << m_numberOfParticles << endl;
    cout << " Number of dimensions : " << m_numberOfDimensions << endl;
    cout << " Number of equilibration Metropolis steps run : 10^" << std::log10(m_equilibrationStepNumber) << endl;
    cout << " Ratio of accepted equilibration steps: " << ((double) m_numberOfAcceptedEquilibrationSteps) / ((double) m_equilibrationStepNumber) << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(m_stepNumber) << endl;
    cout << " Ratio of accepted steps: " << ((double) m_numberOfAcceptedSteps) / ((double) m_stepNumber) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (unsigned int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance in the energy : " << m_energy2 - m_energy * m_energy << endl;
    cout << endl;
}

void Sampler::printOutputToTerminalShort() {
    auto pa = m_waveFunctionParameters;
    auto p = pa.size();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Ratio of accepted equilibration steps: " << ((double) m_numberOfAcceptedEquilibrationSteps) / ((double) m_equilibrationStepNumber) << endl;
    cout << " Ratio of accepted steps: " << ((double) m_numberOfAcceptedSteps) / ((double) m_stepNumber) << endl;
    cout << endl;
    cout << " Number of parameters : " << p << endl;
    for (unsigned int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance in the energy : " << m_energy2 - m_energy * m_energy << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities.
     */
    m_energy = m_cumulativeEnergy / m_stepNumber;
    m_energy2 = m_cumulativeEnergy2 / m_stepNumber;

    //add stuff here to update other sampled quantities for gradient
}


void Sampler::initiateFile(std::string filename){
    std::ofstream file (filename, std::ofstream::app);
    file << "#n_particles" << '\t'
        << "n_dimensions" << '\t'
        << "n_steps" << '\t'
        << "n_accepted_steps" << '\t'
        << "E" << '\t'
        << "var" << '\t'
        << "params" << endl;
    file.close();
}

void Sampler::writeToFile(std::string filename){
    std::ofstream file (filename, std::ofstream::app);
    file << m_numberOfParticles << '\t'
        << m_numberOfDimensions << '\t'
        << m_stepNumber << '\t'
        << m_numberOfAcceptedSteps << '\t'
        << m_energy << '\t'
        << m_energy2 - m_energy * m_energy; // variance
    for (unsigned int i = 0; i < m_waveFunctionParameters.size(); i++)
        file << '\t' << m_waveFunctionParameters[i];
    file << endl;
    file.close();
}