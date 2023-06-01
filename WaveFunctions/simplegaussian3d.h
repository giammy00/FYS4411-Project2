#ifndef __SIMPLE_GAUSSIAN_3D__
#define __SIMPLE_GAUSSIAN_3D__
#include <memory>

#include "wavefunction.h"

class SimpleGaussian3D : public WaveFunction {
public:
    SimpleGaussian3D(double alpha, double beta);
    
    void InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles);
    void adjustPosition(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double> step);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index);
    std::vector<double> quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
    double phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
    std::vector<double> getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles);
};
#endif