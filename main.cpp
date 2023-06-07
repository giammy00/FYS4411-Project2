#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <chrono>
#include <string>
#include <nlopt.hpp>
#include <Math/random.h>
#include "utils.h"
#include "rbmparams.hpp"
using namespace std;
/*
h-bar = 1
m = 1
omega = 1
*/
int main(int argc, char *argv[]) {

    ///////// Variable declarations and default values //////////////////

    //set a maximum number of iterations for gd
    unsigned int  nMaxIter = 3E0;
    //set tolerance for convergence of gd
    double  energyTol = 1E-7;
    //parameters of the simulation
    SimulationParams simPar ;
    //hyperparameters for gradient descent:
    GdParams gd_parameters ;

    if (argc>1){
        // Input filename from command line
        std::string input_filename = argv[1];
        // create filestream to read in inputs from file
        std::ifstream infile;
        infile.open(input_filename);

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
            {gd_parameters.learning_rate = std::stod(value);}
            else if (name == "momentum")
            {gd_parameters.momentum = std::stod(value);}
            else if (name == "gamma")
            {simPar.gamma = std::stod(value);}
            else if (name == "nMaxIter")
            {nMaxIter = std::stoi(value);}
            else if (name == "energyTol")
            {energyTol = std::stod(value);}
            else if (name == "calculateGradients")
            {simPar.calculateGradients = (bool)std::stoi(value);}
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
            else if (name == "sigma")
            {simPar.sigma = std::stod(value);}
            else if (name == "filename")
            {simPar.filename = value;}
            else if (name == "numberOfDimensions")
            {simPar.numberOfDimensions=std::stoi(value);}
            else if (name == "numberOfHiddenNodes")
            {simPar.N_hidden=std::stoi(value);}
            else
            {
                std::cout << "Error reading file." << std::endl;
                return 1;
            }
        }
    }
    else    {
    std::cout << "error: No settings file provided.\n" << 
    "Please use as: " << argv[0] << " <settings_file> [ , <algorithm>].\n" << 
    "Program will be terminated." << std::endl;
    exit(1);
    }
    std::string OPTIMIZATION_ALGORITHM = "LBFGS";
    if (argc>2){
        std::string algo = argv[2];
        if (algo == "GD"){
            OPTIMIZATION_ALGORITHM=algo;
        }
        else if (algo == "LBFGS") {
            ; //do nothing already set, just to display error if not properly set.
        }
        else {
            cout << "Error. Algorithm not recognized.\n Supported algorithms:  'GD' or 'LBFGS', defaults to LBFGS.\n " <<
                "Please use as: " << argv[0] << " <settings_file> [ , <algorithm>].\n" <<  endl;
            exit(1);
        }
    }
    

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

    //
    
    ///// THE FOLLOWING RUNS GRADIENT DESCENT WITH THE TWO IMPLEMENTED METHODS!
    if (simPar.calculateGradients){
        cout << "Optimizing wave function with gradient method " << OPTIMIZATION_ALGORITHM << endl;
        double optimal_energy;
        //this rng is used to construct the wavefunction parameters
        auto rngP =  Random(32384);
        RBMParams parameters = RBMParams(2*simPar.numberOfParticles, //number of visible nodes
                                        simPar.N_hidden,  //number of hidden nodes   
                                        &rngP);
        simPar.rbmParamsPtr=&parameters;
        if (OPTIMIZATION_ALGORITHM=="GD"){
            momentumOptimizer opt(simPar.rbmParamsPtr->m_allParams.size(), &gd_parameters);
            opt.set_min_objective(wrapSimulation, (void *) & simPar );
            opt.set_maxeval(nMaxIter);
            opt.set_ftol_abs(energyTol);
            opt.optimize(parameters.m_allParams, optimal_energy);
        }
        else{
            // LBFGS ALGORITHM BY DEFAULT: 
            nlopt::opt opt(nlopt::LD_LBFGS, simPar.rbmParamsPtr->m_allParams.size());
            opt.set_min_objective(wrapSimulation, (void *) & simPar );
            opt.set_maxeval(nMaxIter);
            opt.set_ftol_abs(energyTol);
            opt.optimize(parameters.m_allParams, optimal_energy);
        }
    }
    //////////////// THE FOLLOWING LINE DOES A LONG SIMULATION WITHOUT GRADIENT COMPUTATION
    else{
        cout << "Running simulation without computing gradients.\n " <<
        "If you wish to enable gradient method, set calculateGradients=0 in simulation_input.txt. " << endl; 
        RBMParams parameters(2*simPar.numberOfParticles, simPar.N_hidden, "./Outputs/optimizedRBMParams.bin");
        simPar.rbmParamsPtr=&parameters;
        wrapSimulationLargeScale( parameters.m_allParams, (void*) & simPar);
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
