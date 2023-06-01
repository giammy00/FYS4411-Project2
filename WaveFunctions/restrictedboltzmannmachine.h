#ifndef __RESTRICTED_BOLTZMANN_MACHINE__
#define __RESTRICTED_BOLTZMANN_MACHINE__
#include<iostream>
#include"wavefunction.h"
#include"../rbmparams.hpp"
#include<memory>
// struct rbmParams{
//     //a struct to contain the parameters of the boltzmann machine
//     std::vector<double> a;
//     std::vector<double> b;
//     std::vector<std::vector<double>> W;
// };
// typedef struct rbmParams RBMParams; 
// RBMParams initWeights(int nVisible, int Nhidden,  Random * rng);

class RestrictedBoltzmannMachine : public WaveFunction {
    /////Implementation of a RBM marginal probability density function
    public:
        RestrictedBoltzmannMachine(double sigma, RBMParams * trainableParameters);//done
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);//done
        void InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles);//done
        void adjustPosition(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double> step) ;//done
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);//done
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index);//done
        std::vector<double> quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);//done
        double phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);//done
        std::vector<double> getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles);//done
        const std::vector<double>& getParameters()   ;
    private:
    //parameters of restricted boltzmann machine
        double m_sigma ;
        double m_sigma2;
        RBMParams * m_trainableParameters;
        unsigned int m_Nhidden, m_Nvisible;
        double m_gaussianTerm;
        double m_productTerm;
        //store proposal for product of terms (1+exp(sj))
        double m_proposalProduct;
        //stores exp(delta s_j)
        std::vector<double> m_proposalProductExpDelta; 
        //store exp(delta gaussian term)
        double m_proposalGaussianExpDelta;
        std::vector<double> m_expBPlusSumXw;//store exp(b_k+m_SumXw)
        std::vector<std::vector<double>> m_inverseDistances;

};
#endif