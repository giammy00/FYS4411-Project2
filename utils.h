#include<vector>
#include <memory>
#include"sampler.h"
//define a struct as a container of simulation params.
struct simulationParams{
    unsigned int numberOfDimensions;
    unsigned int numberOfParticles;
    unsigned int numberOfMetropolisSteps;
    unsigned int numberOfEquilibrationSteps;
    bool calculateGradients;
    double omega;
    double gamma;
    double a_ho;
    double stepLength;
    std::string filename;
};

typedef struct simulationParams SimulationParams;
double wrapSimulation(  
    const std::vector<double> &params,
    std::vector<double> &grad, 
    SimulationParams P  );
std::unique_ptr<Sampler> runSimulation(
    SimulationParams P,
    std::vector<double> params
);