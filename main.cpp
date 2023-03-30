#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <chrono>
#include <string>
#include <nlopt.hpp>
#include "utils.h"
using namespace std;
/*
h-bar = 1
m = 1
omega = 1
*/


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
    unsigned int iterCount, nMaxIter = 3E0;
    //set tolerance for convergence of gd
    double energyChange, oldEnergy, newEnergy, energyTol = 1E-7;
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
