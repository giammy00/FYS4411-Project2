
#include <memory>
#include <vector>
#include <string>
#include <math.h>
#include "restrictedboltzmannmachine.h"
#include "../Math/random.h"

RestrictedBoltzmannMachine::RestrictedBoltzmannMachine(double sigma, int Nvisible, int Nhidden, Random* rng){
//initialize a Restricted Boltzmann machine wavefunction:
//all weights are initialized randomly using gaussian distribution
//Nvisible is M in the project (=P*D)
//Nhidden is N in the project
// THIS IS TENTATIVE . CONSTRUCTOR BETTER NEEDS THE PARAMETER OF THE RBM AS INPUT. (THEY NEED TO BE UPDATED OUTSIDE AND
// PASSED INTO THIS FUNCTION ITERATIVELY TO DO GD)

    m_sigma = sigma;
    m_a.reserve(Nvisible);
    m_b.reserve(Nhidden);
    m_W.reserve(Nvisible);
    double tmp;
    //these factors are used to scale the initial random parameters. Arbitrary choice.
    double norm_a = 1.0/Nvisible;
    double norm_b = 1.0/Nhidden;
    double norm_w = norm_a*norm_b;
    //init weights a
    for(size_t i = 0; i<Nvisible; i++){
        tmp = rng->nextGaussian(0, norm_a);
        m_a.push_back(tmp);
    }
    //init weights b
    for(size_t i = 0; i<Nhidden; i++){
        tmp = rng->nextGaussian(0, norm_b);
        m_b.push_back(tmp);        
    }
    //init weights W
    for(size_t i = 0; i<Nvisible; i++){
        std::vector<double> tmp_arr(Nhidden);
        //fill tmp array with random numbers
        for(size_t j =0; j<Nhidden; j++){
            tmp = rng->nextGaussian(0, norm_w);
            tmp_arr[j] = tmp;
        }
        m_W.push_back(tmp_arr);
    }
}
double RestrictedBoltzmannMachine::evaluate(std::vector<std::unique_ptr<class Particle>>& particles){
    //return Psi . NOT Psi^2

}
