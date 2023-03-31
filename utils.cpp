#pragma once
#include<vector>
#include <memory>

#include "utils.h"
#include "omp.h"
#include"sampler.h"
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


double wrapSimulation(const std::vector<double> &params, std::vector<double> &grad, void * xPtr) {
    //wraps the runSimulation function so that it can be used with nlopt library  
    //check if output file already exists, if so, 
    //set file_initiated to true (avoid printing header line over and over in outputs.txt)
    //just to handle the .txt easier to handle.

    //convert P to struct pointer:
    SimulationParams* P = (SimulationParams *) xPtr;
    bool file_initiated;
    if( FILE * fptr  = fopen(P->filename.c_str(),"r") ){
        fclose(fptr);
        file_initiated = true;
    }
    else
        file_initiated = false;
    int NUM_THREADS = omp_get_max_threads();
    std::cout << "Using " << NUM_THREADS << " threads." << std::endl;
    std::vector< std::unique_ptr< class Sampler >> samplers(NUM_THREADS);

    ///////////////////////////////////////////////
    ///// START PARALLEL REGION //////////////////
    ///////////////////////////////////////////////
    
    #pragma omp parallel
    {
        int thread_number = omp_get_thread_num();
        auto sampler = runSimulation(
                P,
                params //params to init wavefunc
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
        collective_sampler->initiateFile(P->filename);
        file_initiated = true;
    }
    collective_sampler->writeToFile(P->filename);            
    //write to terminal:
    collective_sampler->printOutputToTerminalShort();

    //compute energy difference
    double energy = collective_sampler->getEnergy();


    //sampler computes gradient 
    grad = collective_sampler->computeGradientEtrial();

    return energy;
}


std::unique_ptr<Sampler> runSimulation(
    SimulationParams *P,
    std::vector<double> params
){
    int seed = 2023;
    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    // auto particles = setupNonOverlappingGaussianInitialState(numberOfDimensions, numberOfParticles, *rng, a_ho);
    auto particles = setupNonOverlappingGaussianInitialState(P->numberOfDimensions, P->numberOfParticles, *rng, P->a_ho);
    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(
            // Construct unique_ptr to Hamiltonian
            // std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<HarmonicOscillator3D>(P->omega, P->gamma),
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
            P->calculateGradients);
    
    // Run steps to equilibrate particles
    auto sampler = system->runEquilibrationSteps(
            P->stepLength,
            P->numberOfEquilibrationSteps);

    // Run the Metropolis algorithm
    sampler = system->runMetropolisSteps(
            std::move(sampler),
            P->stepLength,
            P->numberOfMetropolisSteps);

    return sampler;
}