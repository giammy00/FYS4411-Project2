#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    double r2; //store sum of squared positions of particles
    SimpleGaussian(double alpha, std::vector<std::unique_ptr<class Particle>>& particles);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index);
    std::vector<double> quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
    double phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
    std::vector<double> SimpleGaussian::getGradientTerms( double localenergy );
    void SimpleGaussian::updateCachedVariables(std::vector<std::unique_ptr<class Particle>>& particles,  std::vector<double>& step);
};
