#include <memory>
#include <math.h>
#include <assert.h>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

void testDoubleDerivative3d(
    std::vector<std::unique_ptr<class Particle>>& particles,
    class WaveFunction& waveFunction,
    double nabla);

void SimpleGaussian::InitialisePositions(std::vector<std::unique_ptr<class Particle>>&){}

void SimpleGaussian::adjustPosition(std::vector<std::unique_ptr<class Particle>>&, int, std::vector<double>){}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    // Returns Phi, not Phi^2
    double r2 = 0;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        for (unsigned int j = 0; j<position.size(); j++)
            r2 += position[j]*position[j];
    }
    return exp(-1*m_parameters[0]*r2);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    //***************WE RETURN d2/dx2(phi)/phi NOT d2/dx2(phi)*********************

    // The second derivative of exp(-alpha x*x) is exp(-alpha x*x)*(4*alpha*alpha*x*x - 2*alpha)
    
    
    // /* Analytical expression
    double r2 = 0;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        for (unsigned int j = 0; j < particles[i]->getNumberOfDimensions(); j++)
            r2 += position[j]*position[j];
    }
    int n = particles.size() * particles[0]->getNumberOfDimensions();
    double nabla2 = 4*m_parameters[0]*m_parameters[0]*r2 - 2*n*m_parameters[0];

    // testDoubleDerivative3d(particles, *this, nabla2);

    return nabla2;
}


std::vector<double> SimpleGaussian::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************

    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    for (unsigned int i=0; i<force.size(); i++){
        force[i] *= -2*m_parameters[0];
    }
    return force;
}

std::vector<double> SimpleGaussian::quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************

    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    for (unsigned int i=0; i<force.size(); i++){
        force[i] += step[i];
        force[i] *= -2*m_parameters[0];
    }
    return force;
}

double SimpleGaussian::phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    // Calculate (phi(new)/phi(old))**2
    auto pos = particles[index]->getPosition();
    double dr2=0;
    for (unsigned int i=0; i<step.size(); i++)
        // dr2 += (2*pos[i]+step[i])*step[i]; // (pos[i]+step[i])*(pos[i]+step[i])-pos[i]*pos[i]
        dr2 += (pos[i]+step[i])*(pos[i]+step[i])-pos[i]*pos[i];
    return exp(-2*m_parameters[0]*dr2);
}

std::vector<double> SimpleGaussian::getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles){
    double r2 = 0;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        for (unsigned int j = 0; j<position.size(); j++)
            r2 += position[j]*position[j];
    }
    return std::vector<double>{-r2};
}