#pragma once
#include <memory>
#include <string>
#include <iostream>

class Sampler {
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        int numberOfWFParams
        );
    // Sampler(std::vector<Sampler> samplers);


    void sample(bool acceptedStep, class System* system);
    void equilibrationSample(bool acceptedStep);
    void transferWaveFunctionParameters(std::vector<double> parameters);
    void printOutputToTerminal();
    void printOutputToTerminalShort();
    void computeAverages();
    double getEnergy() { return m_energy; }

    void initiateFile(std::string filename);
    void writeToFile(std::string filename);

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_equilibrationStepNumber = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;
    unsigned int m_numberOfAcceptedEquilibrationSteps = 0;
    double m_energy = 0;
    double m_cumulativeEnergy = 0;
    double m_energy2 = 0;
    double m_cumulativeEnergy2 = 0;
    //sampled quantities for computing gradient.
    //for the h.o. we only need one double, might need to be changed to std::vector<double> if more than one quantity needs sampling.
    std::vector<double> m_cumulativeGradientTerms;
    std::vector<double> m_gradientTerms;
    std::vector<double> m_waveFunctionParameters;

};
