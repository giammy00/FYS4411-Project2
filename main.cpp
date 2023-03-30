#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <math.h>
#include <chrono>
#include <string>
#include <nlopt.hpp>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/simplegaussian3d.h"
#include "WaveFunctions/interactinggaussian.h"
#include "WaveFunctions/interactinggaussian3d.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillator3d.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Solvers/metropolisHastings.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"
using namespace std;
#include"omp.h"
/*
h-bar = 1
m = 1
omega = 1
*/
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
    string filename;
};
typedef struct simulationParams SimulationParams;

double wrapSimulation(const std::vector<double> &params, std::vector<double> &grad, SimulationParams P) {
    //wraps the runSimulation function so that it can be used with nlopt library  
    //check if output file already exists, if so, 
    //set file_initiated to true (avoid printing header line over and over in outputs.txt)
    //just to handle the .txt easier to handle.
    bool file_initiated;
    if( FILE * fptr  = fopen(P.filename.c_str(),"r") ){
        fclose(fptr);
        file_initiated = true;
    }
    else
        file_initiated = false;
    int NUM_THREADS = omp_get_max_threads();
    cout << "Using " << NUM_THREADS << " threads." << endl;
    std::vector< std::unique_ptr< class Sampler >> samplers(NUM_THREADS);

    ///////////////////////////////////////////////
    ///// START PARALLEL REGION //////////////////
    ///////////////////////////////////////////////
    
    #pragma omp parallel
    {
        int thread_number = omp_get_thread_num();
        auto sampler = runSimulation(
                P,
                params //params to init wavefunction
                );
        samplers[thread_number]= std::move(sampler);
    }
    ///////////////////////////////////////////////
    ///// END PARALLEL REGION //////////////////
    ///////////////////////////////////////////////


    //gather all simulation results in one sampler
    std::unique_ptr< class Sampler > collective_sampler = std::make_unique< class Sampler >( samplers );
    //write to file
    if(!file_initiated){
        collective_sampler->initiateFile(P.filename);
        file_initiated = true;
    }
    collective_sampler->writeToFile(P.filename);            
    //write to terminal:
    collective_sampler->printOutputToTerminalShort();

    //compute energy difference
    double energy = collective_sampler->getEnergy();


    //sampler computes gradient 
    grad = collective_sampler->computeGradientEtrial();

    return energy
}

std::unique_ptr<Sampler> runSimulation(
    SimulationParams P,
    std::vector<double> params
){
    int seed = 2023;
    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    // auto particles = setupNonOverlappingGaussianInitialState(numberOfDimensions, numberOfParticles, *rng, a_ho);
    auto particles = setupNonOverlappingGaussianInitialState(P.numberOfDimensions, P.numberOfParticles, *rng, P.a_ho);
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(
            // Construct unique_ptr to Hamiltonian
            // std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<HarmonicOscillator3D>(P.omega, P.gamma),
            // Construct unique_ptr to wave function
            // std::make_unique<SimpleGaussian>(params[0]),
            // std::make_unique<SimpleGaussian3D>(params[0], params[1]),
            // std::make_unique<InteractingGaussian>(params[0]),
            std::make_unique<InteractingGaussian3D>(params[0], params[1]),
            // Construct unique_ptr to solver, and move rng
            std::make_unique<MetropolisHastings>(std::move(rng)),
            // std::make_unique<Metropolis>(std::move(rng)),
            // Move the vector of particles to system
            std::move(particles),
            P.calculateGradients);
    
    // Run steps to equilibrate particles
    auto sampler = system->runEquilibrationSteps(
            P.stepLength,
            P.numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    sampler = system->runMetropolisSteps(
            std::move(sampler),
            P.stepLength,
            P.numberOfMetropolisSteps);

    return sampler;
}
int main(int argc, char *argv[]) {

    ///////// Variable declarations and default values //////////////////
    //hyperparameters for gradient descent:
    double learning_rate = 3E-3;
    double momentum = 0.6;
    //store initial trainable parameters of the wave function
    //wfParams is reset to wfParams0 every time a new gradient descent is started
    std::vector<double> wfParams = std::vector<double>{0.5, 1.};
    //define gradient placeholder, to which the gradient will be written
    std::vector<double> gradient = std::vector<double>(wfParams.size());
    //for momentum GD:
    std::vector<double> velocity = std::vector<double>(wfParams.size(), 0.0);
    //set a maximum number of iterations for gd
    unsigned int iterCount, nMaxIter = 1E0;
    //set tolerance for convergence of gd
    double energyChange, oldEnergy, newEnergy, energyTol = 1E-7;
    bool calculateGradients = true;

    SimulationParams simPar ;
    if (argc > 1){
        // Input filename from command line
        std::string input_filename = argv[1];
        std::string path_input = "./Input/";

        // create filestream to read in inputs from file
        std::ifstream infile;
        infile.open(path_input + input_filename + ".txt");

        if (!infile.is_open())
        {
            std::cout << "Error opening file" << std::endl;
            return 1;
        }

        // read in inputs from file

        std::string line;

        while (std::getline(infile, line))
        {
            // loop through lines and sett corect variables. name and value separated by "="
            std::string name = line.substr(0, line.find("="));
            std::string value = line.substr(line.find("=") + 1);
            if (value == "")
            {}
            else if (name == "learning_rate")
            {learning_rate = std::stod(value);}
            else if (name == "momentum")
            {momentum = std::stod(value);}
            else if (name == "alpha")
            {wfParams[0] = std::stod(value);}
            else if (name == "beta")
            {wfParams[1] = std::stod(value);}
            else if (name == "gamma")
            {simPar.gamma = std::stod(value);}
            else if (name == "nMaxIter")
            {nMaxIter = std::stoi(value);}
            else if (name == "energyTol")
            {energyTol = std::stod(value);}
            else if (name == "calculateGradients")
            {calculateGradients = (bool)std::stoi(value);}
            else if (name == "numberOfParticles")
            {simPar.numberOfParticles = std::stoi(value);}
            else if (name == "numberOfMetropolisSteps")
            {simPar.numberOfMetropolisSteps = std::stoi(value);}
            else if (name == "numberOfEquilibrationSteps")
            {simPar.numberOfEquilibrationSteps = std::stoi(value);}
            else if (name == "omega")
            {simPar.omega = std::stod(value);}
            else if (name == "stepLength")
            {simPar.stepLength = std::stod(value);}
            else if (name == "filename")
            {simPar.filename = "Outputs/" + value;}
            else
            {
                std::cout << "Error reading file" << std::endl;
                return 1;
            }
        }
    }
    else{
        std::cout << "WARNING: No settings file provided" << std::endl;
    }

    //DEFINE SIMULATION PARAMETERS
    //NO GRADIENT DESCENT PARAMETERS HERE,
    //NOT EVEN THE TRAINABLE PARAMS of the WF.

    simPar.numberOfDimensions=3;
    simPar.numberOfParticles=10;//50; 100
    simPar.numberOfMetropolisSteps=3E4;
    simPar.numberOfEquilibrationSteps=1E3;
    simPar.calculateGradients=true;
    simPar.omega=1;
    simPar.gamma=1;
    simPar.stepLength=5E-1;
    simPar.filename="Outputs/output.txt";
    simPar.a_ho = std::sqrt(1./simPar.omega); // Characteristic size of the Harmonic Oscillator
    simPar.stepLength *= simPar.a_ho; // Scale the steplength in case of changed omega            

    //#define TIMEING // Comment out turn off timing
    #ifdef TIMEING
    auto times = vector<int>();
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();
    #endif

    iterCount=0;
    energyChange=1; //set to 1 just to enter while loop, should be >= energyTol
    oldEnergy = 1E7;//to enter while loop twice
    while(  (iterCount<nMaxIter) && (energyChange>=energyTol)   )
    {
        newEnergy = wrapSimulation(wfParams, gradient, simPar);
        energyChange = fabs(oldEnergy-newEnergy);
        oldEnergy = newEnergy;
        //update parameters using momentum gd
        for(unsigned int i=0; i<gradient.size(); i++){
            velocity[i] = momentum *  velocity[i] - learning_rate * gradient[i] ;
            wfParams[i] += velocity[i];
        }
        iterCount++;
    }
    #ifdef TIMEING
    auto t2 = high_resolution_clock::now();
    /* Getting number of milliseconds as an integer. */
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    times.push_back(ms_int.count());
    #endif
    #ifdef TIMEING
    cout << "times : " << endl;
    for(unsigned int i = 0; i<times.size(); i++)
        cout << times[i] << endl;
    #endif

    return 0;
}
