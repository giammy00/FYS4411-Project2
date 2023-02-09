#pragma once
#include <memory>

class Sampler {
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength,
        unsigned int numberOfMetropolisSteps);


    void sample(bool acceptedStep, class System* system);
    void printOutputToTerminal(class System& system, unsigned int& equilibrationSteps, unsigned int& acceptedEquilibrationSteps);
    void printOutputToTerminalShort(class System& system, unsigned int& equilibrationSteps, unsigned int& acceptedEquilibrationSteps);
    void computeAverages();
    double getEnergy() { return m_energy; }

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;
    double m_energy = 0;
    double m_cumulativeEnergy = 0;
    double m_energy2 = 0;
    double m_cumulativeEnergy2 = 0;
};
