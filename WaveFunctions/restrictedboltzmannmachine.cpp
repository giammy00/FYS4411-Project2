#include <memory>
#include <vector>
#include <string>
#include <math.h>
#include "restrictedboltzmannmachine.h"
#include "../Math/random.h"
RBMParams initWeights(int Nvisible, int Nhidden,  Random * rng){
    //random initializer for trainable parameters of the restricted boltzmann machine

    RBMParams trainableParams ;
    trainableParams.a.reserve(Nvisible);
    trainableParams.b.reserve(Nhidden);
    trainableParams.W.reserve(Nvisible);
    double tmp;
    //these factors are used to scale the initial random parameters. Arbitrary choice.
    double norm_a = 1.0/Nvisible;
    double norm_b = 1.0/Nhidden;
    double norm_w = norm_a*norm_b;
    //init weights a
    for(size_t i = 0; i<Nvisible; i++){
        tmp = rng->nextGaussian(0, norm_a);
        trainableParams.a.push_back(tmp);
    }
    //init weights b
    for(size_t i = 0; i<Nhidden; i++){
        tmp = rng->nextGaussian(0, norm_b);
        trainableParams.b.push_back(tmp);        
    }
    //init weights W
    for(size_t i = 0; i<Nvisible; i++){
        std::vector<double> tmp_arr(Nhidden);
        //fill tmp array with random numbers
        for(size_t j =0; j<Nhidden; j++){
            tmp = rng->nextGaussian(0, norm_w);
            tmp_arr[j] = tmp;
        }
        trainableParams.W.push_back(tmp_arr);
    }
}

RestrictedBoltzmannMachine::RestrictedBoltzmannMachine(double sigma, RBMParams * trainableParameters){
/*constructor of a Restricted Boltzmann machine wavefunction.
*/
    m_sigma = sigma;
    m_sigma2 = sigma*sigma;//always need sigma2
    m_trainableParameters = trainableParameters;
    m_Nhidden = trainableParameters->b.size();
    m_Nvisible = trainableParameters->a.size();
}
double RestrictedBoltzmannMachine::evaluate(std::vector<std::unique_ptr<class Particle>>& particles){
    //return Psi . NOT Psi^2
    /*
    m_gaussian = exp(-sum( (xi-ai)/2sigma^2 ))
    m_productTerm = the other \Prod term appearing in the wf expression
    They are properly init through call to wavefunction->initialisePositions()
    They are updated at every iteration using waveFunction->adjustPosition in call to solver.step()
    */
    return m_gaussianTerm * m_productTerm ; 
}
void RestrictedBoltzmannMachine::InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles){
    /*
    This function does not really initialise ``positions", rather it initializes some of the cached
    quantities that are useful in the computation of Psi (and of the gradients). 
    */
    m_SumSquares = 0.0;
    m_SumXw = std::vector<double>(m_Nvisible, 0.0);
    m_expSumXw = std::vector<double>(m_Nvisible);
    m_productTerm = 1.0;
    double tmp1 ;
    double X_i ; 
    for(unsigned int i =0 ; i<m_Nvisible; i++){
        //visible node i corresponds to 
        //position i%2 of particle i/2 
        X_i = particles[i%2]->positions[i/2];
        tmp1 = X_i -(m_trainableParameters->a[i]); 
        tmp1*=tmp1;
        m_SumSquares+=tmp1;
        for(unsigned_int j=0; j<m_Nhidden; j++){
            m_SumXw[j]+=X_i * m_trainableParameters->W[i][j];
        }
    }
    for(unsigned_int j=0; j<m_Nhidden; j++){
        m_SumXw[j]/=m_sigma2;
        m_expSumXw[j] = exp(m_trainableParameters->b[j] + m_SumXw[j]);
        m_productTerm*=(1.0+m_expSumXw[j]);
    }
    m_SumSquares/=(2*m_sigma2);
    m_gaussianTerm = exp(m_SumSquares);
}

void RestrictedBoltzmannMachine::adjustPosition(std::vector<std::unique_ptr<class Particle>>& particles, 
                                                int index, std::vector<double> step){

    //need to update  all cached stuff!!!
}



