//optimzier for gradient descent.

//has a:
//pointer to wf (which serves to update the parameters of the wf)
//velocity term (to do momentum gd)

#pragma once
#include<vector>
#include"WaveFunctions/wavefunction.h"

class MomentumGD {
    private:
    WaveFunction * m_waveFuncPtr;
    std::vector<double> m_vel;
    double m_momentum;
    double m_lr;

    public: 
    MomentumGD( double lr, double momentum , WaveFunction * wfPtr );
    void step( std::vector<double> gradient );
}