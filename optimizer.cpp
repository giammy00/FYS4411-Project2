#include"optimizer.h"

MomentumGD::MomentumGD(double lr, double momentum, WaveFunction * wavefuncPtr){
    m_lr = lr;
    m_momentum = momentum; 
    m_waveFuncPtr = wavefuncPtr;
}

void MomentumGD::step( std::vector<double> gradient ){
    m_vel = m_momentum*m_vel + m_lr * gradient ;
}
