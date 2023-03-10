#pragma once

#include <memory>

#include "wavefunction.h"

class InteractingGaussian : public WaveFunction {
public:
    InteractingGaussian(double alpha, double beta = 1, double a = 0.0043);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index);
    std::vector<double> quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
    double phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
    double uPrime(double r);
    double uDoublePrime(double r);
};
