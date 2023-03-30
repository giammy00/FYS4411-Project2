#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <math.h>
#include <chrono>
#include <string>

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
std::unique_ptr<Sampler> runSimulation(
    unsigned int numberOfDimensions,
    unsigned int numberOfParticles,
    unsigned int numberOfMetropolisSteps,
    unsigned int numberOfEquilibrationSteps,
    double omega,
    double gamma,
    double a_ho,
    std::vector<double> params,
    double stepLength
){
    int seed = 2023;
    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    // auto particles = setupNonOverlappingGaussianInitialState(numberOfDimensions, numberOfParticles, *rng, a_ho);
    auto particles = setupNonOverlappingGaussianInitialState(numberOfDimensions, numberOfParticles, *rng, a_ho);
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(
            // Construct unique_ptr to Hamiltonian
            // std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<HarmonicOscillator3D>(omega, gamma),
            // Construct unique_ptr to wave function
            // std::make_unique<SimpleGaussian>(params[0]),
            // std::make_unique<SimpleGaussian3D>(params[0], params[1]),
            // std::make_unique<InteractingGaussian>(params[0]),
            std::make_unique<InteractingGaussian3D>(params[0], params[1]),
            // Construct unique_ptr to solver, and move rng
            std::make_unique<MetropolisHastings>(std::move(rng)),
            // std::make_unique<Metropolis>(std::move(rng)),
            // Move the vector of particles to system
            std::move(particles));
    
    // Run steps to equilibrate particles
    auto sampler = system->runEquilibrationSteps(
            stepLength,
            numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    sampler = system->runMetropolisSteps(
            std::move(sampler),
            stepLength,
            numberOfMetropolisSteps);

    return sampler;
}
int main() {
    // Seed for the random number generator
    // int seed = 2023;
    
    //hyperparameters for gradient descent:
    double learning_rate = 3E-3;
    double momentum = 0.6;
    //store initial trainable parameters of the wave function
    //wfParams is reset to wfParams0 every time a new gradient descent is started
    std::vector<double> wfParams, wfParams0 = std::vector<double>{0.50, 2.86};
    int nParams = wfParams0.size();
    //for momentum GD:
    std::vector<double> velocity = std::vector<double>(nParams, 0.0);
    //set a maximum number of iterations for gd
    unsigned int iterCount, nMaxIter = 1E0;
    //set tolerance for convergence of gd
    double energyChange, oldEnergy, newEnergy, energyTol = 1E-10;

    double numberOfParticles;
    auto numberOfParticlesArray=std::vector<unsigned int>{10};//{10,100,500};
    unsigned int numberOfMetropolisSteps = (unsigned int) 3E4;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1E3;
    double omega = 1.0; // Oscillator frequency.
    double gamma = 2.8284; // Harmonic Oscillator flatness.
    double a_ho = std::sqrt(1./omega); // Characteristic size of the Harmonic Oscillator
    double stepLength = 5E-1; // Metropolis step length.
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
                            

    // #define TIMEING // Comment out turn off timing
    #ifdef TIMEING
    auto times = vector<int>();
    #endif

    
    for (unsigned int i = 0; i < numberOfParticlesArray.size(); i++){
        numberOfParticles = numberOfParticlesArray[i];
        iterCount=0;
        energyChange=1; //set to 1 just to enter while loop, should be >= energyTol
        oldEnergy = 1E7;//to enter while loop twice
        wfParams=wfParams0;//restart gradient descent.

        int NUM_THREADS = omp_get_max_threads();
        cout << "Using " << NUM_THREADS << " threads." << endl;
        std::vector< std::unique_ptr< class Sampler > > samplers(NUM_THREADS);

        while(  (iterCount<nMaxIter) && (energyChange>=energyTol)   ){
            
            #ifdef TIMEING
            using std::chrono::high_resolution_clock;
            using std::chrono::duration_cast;
            using std::chrono::duration;
            using std::chrono::milliseconds;
            auto t1 = high_resolution_clock::now();
            #endif
            
            ///////////////////////////////////////////////
            ///// START PARALLEL REGION //////////////////
            ///////////////////////////////////////////////
            
            #pragma omp parallel
            {
                int thread_number = omp_get_thread_num();
                auto sampler = runSimulation(
                        3,//dimensions of our simulation
                        numberOfParticles,
                        numberOfMetropolisSteps,
                        numberOfEquilibrationSteps,
                        omega,
                        gamma,
                        a_ho,
                        wfParams, //I assumed that wave function will take a std::vector<double> of params to be initialized
                        stepLength);
                samplers[thread_number]= std::move(sampler);
            }
            ///////////////////////////////////////////////
            ///// END PARALLEL REGION //////////////////
            ///////////////////////////////////////////////

            #ifdef TIMEING
            auto t2 = high_resolution_clock::now();
            /* Getting number of milliseconds as an integer. */
            auto ms_int = duration_cast<milliseconds>(t2 - t1);
            times.push_back(ms_int.count());
            #endif

            //gather all simulation results in one sampler
            std::unique_ptr< class Sampler > collective_sampler = std::make_unique< class Sampler >( samplers );
            //write to file
            if(!file_initiated){
                collective_sampler->initiateFile(filename);
                file_initiated = true;
            }
            collective_sampler->writeToFile(filename);            
            //write to terminal:
            collective_sampler->printOutputToTerminalShort();

            //compute energy difference
            newEnergy = collective_sampler->getEnergy();
            energyChange = fabs(oldEnergy-newEnergy);
            oldEnergy = newEnergy;

            //sampler computes gradient 
            std::vector<double> gradient = collective_sampler->computeGradientEtrial();

            //update parameters using momentum gd
            for(unsigned int i=0; i<gradient.size(); i++){
                velocity[i] = momentum *  velocity[i] - learning_rate * gradient[i] ;
                wfParams[i] += velocity[i];
            }
            iterCount++;
        }
    }
    #ifdef TIMEING
    cout << "times : " << endl;
    for(unsigned int i = 0; i<times.size(); i++)
        cout << times[i] << endl;
    #endif

    return 0;
}
