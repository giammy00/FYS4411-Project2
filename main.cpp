#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <math.h>
#include <chrono>
#include <string>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Solvers/metropolisHastings.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"
using namespace std;

/*
h-bar = 1
m = 1
omega = 1
*/
std::unique_ptr<Sampler> runSimulation(
    unsigned int numberOfDimensions,
    unsigned int numberOfParticles,
    unsigned int numberOfMetropolisSteps,
    unsigned int numberOfEquilibrationSteps,
    double omega,
    double a_ho,
    std::vector<double> params,
    double stepLength
){
    int seed = 2023;
    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomGaussianInitialState(numberOfDimensions, numberOfParticles, *rng, a_ho);//[x]
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(//[x]
            // Construct unique_ptr to Hamiltonian
            std::make_unique<HarmonicOscillator>(omega),//[x]
            // Construct unique_ptr to wave function
            std::make_unique<SimpleGaussian>(params[0], particles),//[x]
            // Construct unique_ptr to solver, and move rng
            std::make_unique<MetropolisHastings>(std::move(rng)),//[x]
            // std::make_unique<Metropolis>(std::move(rng)),//[x]
            // Move the vector of particles to system
            std::move(particles));
    
    // Run steps to equilibrate particles
    auto sampler = system->runEquilibrationSteps(//[x]
            stepLength,
            numberOfEquilibrationSteps);
    
    // Run the Metropolis algorithm
    sampler = system->runMetropolisSteps(//[ ]
            std::move(sampler),
            stepLength,
            numberOfMetropolisSteps);
    
    //here should call sampler->computeGradientEtrial()
    // here also print info about current energy and variational parameter (move from main)
    // then compute optimizer->step()


    return sampler;
}
int main() {
    // Seed for the random number generator
    // int seed = 2023;
    
    //hyperparameters for gradient descent:
    double learning_rate = 0.001;
    double momentum = 0.9;
    //set initial trainable parameters of the wave function
    std::vector<double> wfParams = std::vector<double>{0.1};
    int nParams = wfParams.size();
    //for momentum GD:
    std::vector<double> velocity = std::vector<double>(nParams, 0.0);
    //set a maximum number of iterations for gd
    unsigned int nMaxIter = 10;
    unsigned int iterCount;
    //set tolerance for convergence of gd
    double energyTol = 0.01;
    double energyChange;
    double oldEnergy, newEnergy; 

    unsigned int numberOfParticles = 1;
    auto numberOfParticlesArray=std::vector<unsigned int>{1,10,100,500};
    unsigned int numberOfMetropolisSteps = (unsigned int) 1E6;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1E5;
    double omega = 1.0; // Oscillator frequency.
    double a_ho = std::sqrt(1./omega); // Characteristic size of the Harmonic Oscillator
    // double alpha = 0.5; // Variational parameter.
    double stepLength = 5E-2; // Metropolis step length.
    stepLength *= a_ho; // Scale the steplength in case of changed omega
    string filename = "Outputs/output.txt";

    //check if file already exists, if so, set file_initiated to true (avoid printing header line over and over in outputs.txt)
    //just to handle the .txt easier to handle.
    bool file_initiated;
    if( FILE * fptr  = fopen(filename.c_str(),"r") ){
        fclose(fptr);
        file_initiated = true;
    }
    else
        file_initiated = false;

    #define TIMEING // Comment out turn off timing
    #ifdef TIMEING
    auto times = vector<int>();
    #endif

    
    for (unsigned int numberOfDimensions = 3; numberOfDimensions < 4; numberOfDimensions++){
        for (unsigned int i = 0; i < numberOfParticlesArray.size(); i++){
            numberOfParticles = numberOfParticlesArray[i];
            iterCount=0;
            energyChange=1; //set to 1 just to enter while loop, should be >= energyTol
            oldEnergy = 1E7;
            while(  (iterCount<nMaxIter) & (energyChange>=energyTol)   ){

                #ifdef TIMEING
                using std::chrono::high_resolution_clock;
                using std::chrono::duration_cast;
                using std::chrono::duration;
                using std::chrono::milliseconds;
                auto t1 = high_resolution_clock::now();
                #endif
                
                auto sampler = runSimulation(
                        numberOfDimensions,
                        numberOfParticles,
                        numberOfMetropolisSteps,
                        numberOfEquilibrationSteps,
                        omega,
                        a_ho,
                        wfParams, //I assume that wave function will take a vector<double> of params to be initialized
                        stepLength);
                

                #ifdef TIMEING
                auto t2 = high_resolution_clock::now();
                /* Getting number of milliseconds as an integer. */
                auto ms_int = duration_cast<milliseconds>(t2 - t1);
                times.push_back(ms_int.count());
                #endif

                if(!file_initiated){
                    sampler->initiateFile(filename);
                    file_initiated = true;
                }
                sampler->writeToFile(filename);
                // Output information from the simulation
                sampler->printOutputToTerminalShort();
                
                newEnergy = sampler->getEnergy();
                energyChange = fabs(oldEnergy-newEnergy);
                oldEnergy = newEnergy;

                //sampler computes gradient 
                std::vector<double> gradient = sampler->computeGradientEtrial();

                //update parameters using momentum gd
                for(int i=0; i<nParams; i++){
                    velocity[i] = momentum *  velocity[i] + learning_rate * gradient[i] ;
                    wfParams[i]   -= velocity[i];
                    std::cout << gradient[i] << std::endl;
                    std::cout << velocity[i] << std::endl;
                }

            }
        }
    }
    #ifdef TIMEING
    cout << "times : " << endl;
    for(unsigned int i = 0; i<times.size(); i++)
        cout << times[i] << endl;
    #endif

    return 0;
}
