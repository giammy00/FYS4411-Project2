#pragma once
#include <memory>

class Sampler {
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength
        );


    void sample(bool acceptedStep, class System* system);
    void equilibrationSample(bool acceptedStep);
    void transferWaveFunctionParameters(std::vector<double> parameters);
    void printOutputToTerminal();
    void printOutputToTerminalShort();
    void computeAverages();
    double getEnergy() { return m_energy; }

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
    std::vector<double> m_waveFunctionParameters;
};
