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
        unsigned int numberOfDimensions,
        int numberOfWFParams
        )
{
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    m_cumulativeGradientTerms = std::vector<std::vector<double>>(numberOfWFParams, std::vector<double>(2, 0.0)) ;
    m_gradientTerms = std::vector<std::vector<double>>(numberOfWFParams, std::vector<double>(2, 0.0)) ;

}

// Sampler::Sampler(std::vector<Sampler> samplers){
//     m_numberOfParticles = samplers[0].m_numberOfParticles;
//     m_numberOfDimensions = samplers[0].m_numberOfDimensions;
//     m_cumulativeGradientTerms = std::vector<std::vector<double>>(numberOfWFParams, std::vector<double>(2, 0.0)) ;

// }

void Sampler::equilibrationSample(bool acceptedStep){
    m_equilibrationStepNumber++;
    m_numberOfAcceptedEquilibrationSteps += acceptedStep;
}

void Sampler::sample(bool acceptedStep, System* system) {
    /*sample all the interesting things 
     */
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += localEnergy * localEnergy;

    //sample quantities for gradient computation
    std::vector<double> currentDerivatives =  system->getdPhi_dParams();
    std::vector<std::vector<double>>currentGradientTerms =  computeGradientTerms( currentDerivatives, localEnergy);

    //now currentGradientTerms is a matrix, different rows correspond to different parameters
    //update each of the rowas adding current sampled quantity
    for(unsigned int i=0; i<currentGradientTerms.size();i++){
        for(unsigned int j=0; j<2;j++){
            m_cumulativeGradientTerms[i][j]+=currentGradientTerms[i][j];
        }
    }
    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;

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

    for(size_t i=0; i<m_cumulativeGradientTerms.size();i++){
        for(size_t j=0; j<2;j++){
            std::cout << i << "  " << j << std::endl;
            m_gradientTerms[i][j]=m_cumulativeGradientTerms[i][j]/m_stepNumber;
        }
    }
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

std::vector<std::vector<double>> Sampler::computeGradientTerms( std::vector<double> dPhi_dParams, double Elocal ){
    ////*****NOTE : by dPhi_dParams we mean the gradient of log(wavefunc.) wrt the variational parameters
    int N = dPhi_dParams.size();
    std::vector<std::vector<double>> gradTerms(N, std::vector<double>(2));

    //in every row of the matrix we have terms related to a different variational parameter
    for(int i = 0 ; i<N; i++){
        gradTerms[i][0] = dPhi_dParams[i];
        gradTerms[i][1] = dPhi_dParams[i]*Elocal;
    }
    return gradTerms;
}

//gradient of the trial energy
std::vector<double> Sampler::computeGradientEtrial()
{
int N = m_gradientTerms.size();
std::vector<double> gradient = std::vector<double>(N); 
    for(int i=0; i<N; i++){
        for(int j=0; j<2; j++){
            gradient[i] = 2*(m_gradientTerms[i][1]-m_energy*m_gradientTerms[i][0]);
        }
    }
return gradient ; 
}