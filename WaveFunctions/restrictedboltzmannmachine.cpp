#include <memory>
#include <vector>
#include <string>
#include <math.h>
#include "restrictedboltzmannmachine.h"
#include "../Math/random.h"
#include "../particle.h"

RBMParams initWeights(int Nvisible, int Nhidden,  Random * rng){
    //random initializer for trainable parameters of the restricted boltzmann machine
    //parameters are stored in a struct 
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
    They are updated at every metropolis step using waveFunction->adjustPosition in call to solver.step()
    */
    return m_gaussianTerm * m_productTerm ; 
}
void RestrictedBoltzmannMachine::InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles){
    /*
    This function does not really initialise ``positions", rather it initializes some of the cached
    quantities that are useful in the computation of Psi (and of the gradients). The cached quantities are:
        - (std::vector) m_expBPlusSumXw  ---> contains at index j:  exp( b_j + \sum_i X_i W_{ij}   ) 
        - double m_productTerm ---> contains \Prod_j (1 + exp( m_expBPlusSum[j]) ) 
        - double m_gaussianTerm  ---> contains exp[ - \sum_i (X_i - a_i)/(2*sigma^2) ]
    */

    double sumSquares = 0.0;
    //will contain b_j + \sum_i X_i W_{ij} 
    std::vector<double> bPlusSumXw = std::vector<double>(m_Nhidden, 0.0)  ;
    m_expBPlusSumXw = std::vector<double>(m_Nhidden);
    m_productTerm = 1.0;
    double tmp1 ;
    double X_i ; 
    for(unsigned int i =0 ; i<m_Nvisible; i++){
        //visible node i corresponds to 
        //position i%2 of particle i/2 
        X_i = particles[i/2]->positions[i%2];
        tmp1 = X_i -(m_trainableParameters->a[i]); 
        tmp1*=tmp1;
        sumSquares+=tmp1;
        for(unsigned int j=0; j<m_Nhidden; j++){
            bPlusSumXw[j]+=X_i * m_trainableParameters->W[i][j];
        }
    }

    for(unsigned int j=0; j<m_Nhidden; j++){
        bPlusSumXw[j]/=m_sigma2;
        bPlusSumXw[j]+=m_trainableParameters->b[j];
        m_expBPlusSumXw[j] = exp( bPlusSumXw[j]);
        m_productTerm*=(1.0+m_expBPlusSumXw[j]);
    }
    sumSquares/=(2*m_sigma2);
    m_gaussianTerm = exp(-sumSquares);
}

void RestrictedBoltzmannMachine::adjustPosition(std::vector<std::unique_ptr<class Particle>>& particles, 
                                                int index, std::vector<double> step)
    //note. Particle index , position[i] (i=0,1) corresponds to visible node of index 2*index+i
{
    //need to update  all cached stuff!!!
    //update product term and cached exp
    m_productTerm = 1.0; 
    //productTerm needs actually to be recomputed from scratch, but useful to compute it here,
    //because there is here a loop over hidden nodes ....could also compute the product term in .evaluate()?
    //yes, but it would take another loop over hidden nodes. so we just exploit the loop which is here.
                        
    double delta ;
    for(unsigned int j=0; j<m_Nhidden; j++){
        //change in the \sum_i X_i w_{ij}
        delta = (   step[0]*m_trainableParameters->W[2*index][j]+
                    step[1]*m_trainableParameters->W[2*index+1][j]  ) /m_sigma2;//delta_k *w_ik /sigma_^2
        m_expBPlusSumXw[j]*=exp(delta);
        m_productTerm*=(1.0+ m_expBPlusSumXw[j]);
    }
    
    //update gaussian term
    std::vector<double> X = particles[index]->getPosition(); 
    delta = step[0]*(  step[0]+2.0*(X[0]-m_trainableParameters->a[2*index])   ) + 
            step[1]*(  step[1]+2.0*(X[1]-m_trainableParameters->a[2*index+1]) );
    m_gaussianTerm*=exp(  - delta )  ;
}

//NOTE: THERE IS A CACHE i.e. some of the quantities needed for the computation are already stored in the class
//(see functions above)
double RestrictedBoltzmannMachine::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles){
    //compute the laplacian of the wave function (see ipynb on boltzmann machines)
    // we use the expression (108) from the lecture notes on boltzmann machines
    double sum=0.0;
    double W_ij;
    for(unsigned int i =0 ; i<m_Nvisible; i++){
        for (unsigned int j=0; j<m_Nhidden; j++){
            W_ij=m_trainableParameters->W[i][j];
            sum+=W_ij*W_ij*(m_expBPlusSumXw[j])/((1.0+m_expBPlusSumXw[j])*(1.0+m_expBPlusSumXw[j]));
        }
        
    }
    sum/=m_sigma2*m_sigma2;
    return -m_Nvisible/m_sigma2+sum;
}
 
std::vector<double> RestrictedBoltzmannMachine::getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles){
    //compute the derivative of log(Psi) wrt variational parameters (see same notes, again)
    //note that the parameters are stored in a struct containing three std::vector! BUT we should return 
    //ONE flattenened vector , like : [ d/da1, d/da2, ..., d/daM, d/db1..., d/dbN, d/dW11, d/dW12...., d/dWMN ]......
    // easiest way: loop through the vectors in the struct , one by one, compute the derivative wrt that parameter
    // append it to the vec with .push_back( )
    std::vector<double> grad;
    std::vector<double> X ;
    double tm1;
    double tm2;
    for (unsigned int i =0 ; i<m_Nvisible/2; i++){
        X = particles[i]->getPosition();
        tm1=X[0]-m_trainableParameters->a[i];
        tm1/=m_sigma2;
        tm2=X[1]-m_trainableParameters->a[i];
        tm2/=m_sigma2;
        grad.push_back(tm1+tm2);
    }
    for (unsigned int i =0 ; i<m_Nhidden; i++){
        tm1=1.0/m_expBPlusSumXw[i]+1.0;
        grad.push_back(1.0/tm1);
    }
    for (unsigned int i =0 ; i<m_Nvisible; i++){
        X = particles[i]->getPosition();
        for (unsigned int j =0 ; j<m_Nhidden; j++){
            tm1=1.0/m_expBPlusSumXw[j]+1.0;
            tm1*=m_sigma2;
            grad.push_back((X[0]+X[1])/tm1);
        }
    }
    return grad;
}



std::vector<double> RestrictedBoltzmannMachine::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index){
    //compute quantum force (see photo on whatsapp)
    std::vector<double> QF;
    std::vector<double> X;
    double a_0, a_1;
    double W_0j, W_1j;
    double sum;
    double delta;
    sum=0.0;
    X= particles[index]->getPosition();
    a_0=m_trainableParameters->a[0+index];
    a_1=m_trainableParameters->a[1+index];
    delta=-(X[0]-a_0)/m_sigma2-(X[1]-a_1)/m_sigma2;
    for(unsigned int j =0 ; j<m_Nhidden; j++){
        W_0j=m_trainableParameters->W[index+0][j];
        W_1j=m_trainableParameters->W[index+1][j];
        sum+=(W_0j+W_1j)/(1.0+1.0/m_expBPlusSumXw[j]);
    }
    QF.push_back(delta+sum);
    return QF;
}

std::vector<double> RestrictedBoltzmannMachine::quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){

    std::vector<double> QF;
    std::vector<double> X;
    double a_0, a_1;
    double W_0j, W_1j;
    double sum;
    double delta;
    sum=0.0;
    X= particles[index]->getPosition();
    a_0=m_trainableParameters->a[0+index];
    a_1=m_trainableParameters->a[1+index];
    delta=-(X[0]+step[0]-a_0)/m_sigma2-(X[1]+step[1]-a_1)/m_sigma2;
    for(unsigned int j =0 ; j<m_Nhidden; j++){
        W_0j=m_trainableParameters->W[index+0][j];
        W_1j=m_trainableParameters->W[index+1][j];
        sum+=W_0j/(1.0+1.0/m_expBPlusSumXw[j]*exp(-1/m_sigma2*W_0j*step[0]))+W_1j/(1.0+1.0/m_expBPlusSumXw[j]*exp(-1/m_sigma2*W_1j*step[1]));
    }
    QF.push_back(delta+sum);
    return QF;
}

double RestrictedBoltzmannMachine::phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    // computated by taking the (phi(new)/phi(old))**2 from the expression (32) in the lecture note
    double A;
    double a_0, a_1;
    std::vector<double> X;
    a_0=m_trainableParameters->a[0+index];
    a_1=m_trainableParameters->a[1+index];
    A=exp((step[0]*(a_0-X[0])+step[1]*(a_1-X[1]))/m_sigma2);

    
}