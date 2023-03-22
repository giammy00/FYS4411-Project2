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
#include "optimizer.h"

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
    double alpha,
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
            std::make_unique<SimpleGaussian>(alpha),//[x]
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
    
    return sampler;
}
int main() {
    // Seed for the random number generator
    // int seed = 2023;
    
    // unsigned int numberOfDimensions = 3;

    //set a maximum number of iterations for gd
    unsigned int nMaxIter = 10;
    unsigned int iterCount;
    //set tolerance for convergence of gd
    double energyTol = 0.01;
    double energyChange;
    double oldEnergy, newEnergy; 
    //set initial parameter for gd
    double alpha;
    //init optimizer for gd
    MomentumGD optimizer = MomentumGD();
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
    for (unsigned int numberOfDimensions = 1; numberOfDimensions < 4; numberOfDimensions++){
        for (unsigned int i = 0; i < numberOfParticlesArray.size(); i++){
            numberOfParticles = numberOfParticlesArray[i];
            //while ( convergence criterion not met )
            iterCount=0;
            alpha=0.3;
            energyChange=1; //set to 1 just to enter while loop, should be >= energyTol
            
            while(iterCount<nMaxIter & energyChange>=energyTol){

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
                        alpha,
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
