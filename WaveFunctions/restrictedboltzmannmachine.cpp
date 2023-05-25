#include <memory>
#include <vector>
#include <string>
#include <math.h>
#include<cassert>
#include "restrictedboltzmannmachine.h"
#include "../Math/random.h"
#include "../particle.h"


RestrictedBoltzmannMachine::RestrictedBoltzmannMachine(double sigma, RBMParams * trainableParameters){
/*constructor of a Restricted Boltzmann machine wavefunction.
*/
    m_sigma = sigma;
    m_sigma2 = sigma*sigma;//always need sigma2
    m_trainableParameters = trainableParameters;
    m_Nhidden = trainableParameters->m_Nhidden;
    m_Nvisible = trainableParameters->m_Nvisible;
    m_numberOfParameters = m_Nvisible+m_Nhidden+m_Nvisible*m_Nhidden;
    m_proposalProductExpDelta = std::vector<double>(m_Nhidden, 0.0);
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
        X_i = particles[i/2]->getPosition()[i%2];
        tmp1 = X_i -(m_trainableParameters->m_a[i]); 
        tmp1*=tmp1;
        sumSquares+=tmp1;
        for(unsigned int j=0; j<m_Nhidden; j++){
            bPlusSumXw[j]+=X_i * m_trainableParameters->m_W[i][j];
        }
    }

    for(unsigned int j=0; j<m_Nhidden; j++){
        bPlusSumXw[j]/=m_sigma2;
        bPlusSumXw[j]+=m_trainableParameters->m_b[j];
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
    //m_productTerm = 1.0; 
    //productTerm needs actually to be recomputed from scratch, but useful to compute it here,
    //because there is here a loop over hidden nodes ....could also compute the product term in .evaluate()?
    //yes, but it would take another loop over hidden nodes. so we just exploit the loop which is here.
                        
    for(unsigned int j=0; j<m_Nhidden; j++){
    //change in the \sum_i X_i w_{ij}
        // delta = (   step[0]*m_trainableParameters->W[2*index][j]+
        //             step[1]*m_trainableParameters->W[2*index+1][j]  ) /m_sigma2;//delta_k *w_ik /sigma_^2
        m_expBPlusSumXw[j]*=m_proposalProductExpDelta[j];//it is exp(delta(s_j) )
    }
    //
    // was already computed in phiratio!
    m_productTerm=m_proposalProduct;
    //update gaussian term
    // std::vector<double> X = particles[index]->getPosition(); 
    // delta = step[0]*(  step[0]+2.0*(X[0]-m_trainableParameters->a[2*index])   ) + 
    //         step[1]*(  step[1]+2.0*(X[1]-m_trainableParameters->a[2*index+1]) );
    m_gaussianTerm*= m_proposalGaussianExpDelta   ; //it is exp(delta(gaussian exponent))

}

//NOTE: THERE IS A CACHE i.e. some of the quantities needed for the computation are already stored in the class
//(see functions above)
double RestrictedBoltzmannMachine::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles){
    //compute the laplacian of Psi (divided by Psi)
    // we use the expression (108) from the lecture notes on boltzmann machines
    // the factor -1/2 is multiplied by the hamiltonian. 
    double secondDerivative;
    double firstDerivative;
    double sum1=0.0; 
    double sum2=0.0; //store sum of 1st and 2nd derivatives
    double W_ij;
    double x_i ;
    for(unsigned int i =0 ; i<m_Nvisible; i++){

        //Compute derivatives of ln(Psi) wrt x_i
        firstDerivative = 0.0;
        secondDerivative = 0.0;
        for (unsigned int j=0; j<m_Nhidden; j++){
            W_ij=m_trainableParameters->m_W[i][j];
            secondDerivative+=W_ij*W_ij*(m_expBPlusSumXw[j])/((1.0+m_expBPlusSumXw[j])*(1.0+m_expBPlusSumXw[j]));
            firstDerivative+= W_ij*m_expBPlusSumXw[j]/(m_expBPlusSumXw[j]+1);
        }
        x_i = particles[i/2]->getPosition()[i%2];
        firstDerivative-=(x_i-m_trainableParameters->m_a[i]);
        sum1+=firstDerivative*firstDerivative;
        sum2+=secondDerivative;
    }
    // NEED ALSO THE FIRST DERIVATIVE OF THE LOG
    sum1/=m_sigma2*m_sigma2;
    sum2/=m_sigma2*m_sigma2;
    sum2+= -static_cast<double>(m_Nvisible)/m_sigma2;
    return sum1+sum2;
}
 
std::vector<double> RestrictedBoltzmannMachine::getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles){
    //compute the derivative of log(Psi) wrt variational parameters (see same notes, again)
    //note that the parameters are stored in a struct containing three std::vector! BUT we should return 
    //ONE flattenened vector , like : [ d/da1, d/da2, ..., d/daM, d/db1..., d/dbN, d/dW11, d/dW12...., d/dWMN ]......
    // easiest way: loop through the vectors in the struct , one by one, compute the derivative wrt that parameter
    // append it to the vec with .push_back( )
    std::vector<double> grad;
    grad.reserve(m_Nhidden+m_Nvisible+m_Nhidden*m_Nvisible);
    std::vector<double> X ;
    double tm1;
    //d/da_i
    for (unsigned int i =0 ; i<m_Nvisible; i++){
        X = particles[i/2]->getPosition();
        tm1=X[i%2]-m_trainableParameters->m_a[i];
        tm1/=m_sigma2;
        grad.push_back(0.5*tm1);
    }
    //d/db_j
    for (unsigned int i =0 ; i<m_Nhidden; i++){
        tm1=1.0/m_expBPlusSumXw[i]+1.0;
        grad.push_back(0.5/tm1);
    }
    //d/W_ij
    for (unsigned int i =0 ; i<m_Nvisible; i++){
        X = particles[i/2]->getPosition();
        for (unsigned int j =0 ; j<m_Nhidden; j++){
            tm1=1.0/m_expBPlusSumXw[j]+1.0;
            tm1*=m_sigma2;
            grad.push_back(0.5*X[i%2]/tm1);
        }
    }
    return grad;
}



std::vector<double> RestrictedBoltzmannMachine::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index){
    //computes quantum force at current position
    std::vector<double> QF = std::vector<double>(2);
    std::vector<double> X;
    double a_0, a_1;
    double W_0j, W_1j;
    double sum1, sum0;
    double delta0, delta1;
    sum1=0.0;
    sum0=0.0;
    X= particles[index]->getPosition();
    a_0=m_trainableParameters->m_a[0+2*index];
    a_1=m_trainableParameters->m_a[1+2*index];
    delta0 =-(X[0]-a_0)/m_sigma2;
    delta1 = -(X[1]-a_1)/m_sigma2;
    for(unsigned int j =0 ; j<m_Nhidden; j++){
        double tmp = (1.0+1.0/m_expBPlusSumXw[j]);
        W_0j=m_trainableParameters->m_W[2*index+0][j];
        W_1j=m_trainableParameters->m_W[2*index+1][j];
        sum0+=(W_0j)/tmp;
        sum1+=(W_1j)/tmp;
    }
    QF[0]=delta0+sum0;
    QF[1]=delta1+sum1;
    return QF;
}

std::vector<double> RestrictedBoltzmannMachine::quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    //computes the quantum force in the proposal position during metropolis hastings.
    std::vector<double> QF;
    QF.reserve(2);
    std::vector<double> X;
    double a_0, a_1;
    double W_0j, W_1j;
    double sum1, sum0;
    double tmp_exp;
    sum0=0.0;
    sum1=0.0;
    X= particles[index]->getPosition();
    a_0=m_trainableParameters->m_a[0+2*index];
    a_1=m_trainableParameters->m_a[1+2*index];
    sum0-=(X[0]+step[0]-a_0) ;
    sum1-=(X[1]+step[1]-a_1) ;

    for(unsigned int j =0 ; j<m_Nhidden; j++){
        W_0j=m_trainableParameters->m_W[2*index+0][j];
        W_1j=m_trainableParameters->m_W[2*index+1][j];
        tmp_exp = m_expBPlusSumXw[j]*m_proposalProductExpDelta[j];
        sum0+=W_0j/(1.0+1.0/tmp_exp);
        sum1+=W_1j/(1.0+1.0/tmp_exp);
    }
    //need a vector of two elements
    QF.push_back(sum0/m_sigma2);
    QF.push_back(sum1/m_sigma2);
    return QF;
}

double RestrictedBoltzmannMachine::phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    // computes (psi(new)/psi(old))**2
    double A;
    double a_0, a_1;
    std::vector<double> X = particles[index]->getPosition();
    a_0=m_trainableParameters->m_a[0+2*index];
    a_1=m_trainableParameters->m_a[1+2*index];
    A=  -step[0]*(2*(X[0]-a_0)+step[0])
        -step[1]*(2*(X[1]-a_1)+step[1]);
    A=exp(A);
    //cache proposal quantities, used in adjustPosition (no need to recompute them!!)
    m_proposalGaussianExpDelta = A;
    m_proposalProduct=1.0;//1+exp(bj+sum xiwij/sigma_2+step*wkj/sigma_2)
    double W_0j, W_1j;
    for (unsigned int j =0 ; j<m_Nhidden; j++){
        W_0j=m_trainableParameters->m_W[2*index+0][j];
        W_1j=m_trainableParameters->m_W[2*index+1][j];
        m_proposalProductExpDelta[j] = exp( (step[0]*W_0j+step[1]*W_1j)/m_sigma2);
        m_proposalProduct*=1.0+m_expBPlusSumXw[j]*m_proposalProductExpDelta[j];
    }
    A*=m_proposalProduct/m_productTerm;
    A*=A;//want the square of psi, not psi
    return A;
}

const std::vector<double>& RestrictedBoltzmannMachine::getParameters() {
    return m_trainableParameters->m_allParams;
}