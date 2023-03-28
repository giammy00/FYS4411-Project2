#include <memory>
#include <math.h>
#include <assert.h>

#include "simplegaussian3d.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

void testDoubleDerivative3d(
    std::vector<std::unique_ptr<class Particle>>& particles,
    class WaveFunction& waveFunction,
    double nabla);

SimpleGaussian3D::SimpleGaussian3D(double alpha, double beta)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
}

void SimpleGaussian3D::InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles){
    assert(particles[0]->getNumberOfDimensions()==3);
}

void SimpleGaussian3D::adjustPosition(std::vector<std::unique_ptr<class Particle>>&, int, std::vector<double>){}

double SimpleGaussian3D::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    // Returns Phi, not Phi^2
    double r2 = 0, beta = m_parameters[1];
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        r2 += position[0]*position[0];
        r2 += position[1]*position[1];
        r2 += beta*position[2]*position[2];
    }
    return exp(-1*m_parameters[0]*r2);
}

double SimpleGaussian3D::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
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
    double r2 = 0, alpha = m_parameters[0], beta = m_parameters[1];
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        r2 += position[0]*position[0];
        r2 += position[1]*position[1];
        r2 += beta*beta*position[2]*position[2]; // Beta squared, because alpha is squared later too
    }
    double nabla2 = 4*alpha*alpha*r2 - (4+2*beta)*particles.size()*alpha; // 4*alpha + 2*alpha*beta

    // TESTING:
    // testDoubleDerivative3d(particles, *this, nabla2);

    return nabla2;
}

// void testDoubleDerivative3d(
//     std::vector<std::unique_ptr<class Particle>>& particles,
//     class WaveFunction& waveFunction,
//     double nablaAnal){
    
//     // Numerical calculation
//     double nabla2 = 0, phi, phi_plus, phi_minus;
//     // double phi, phi_plus, phi_minus;
//     // nabla2 = 0;
//     const double dx = 1e-5, dx2_1 = 1/(dx*dx);
//     phi = waveFunction.evaluate(particles);
//     for (unsigned int i = 0; i < particles.size(); i++){
//         for (unsigned int j = 0; j < 3; j++){
//             particles[i]->adjustPosition(dx, j);
//             phi_plus = waveFunction.evaluate(particles);
//             particles[i]->adjustPosition(-2*dx, j);
//             phi_minus = waveFunction.evaluate(particles);
//             particles[i]->adjustPosition(dx, j);
            
//             nabla2 += (phi_plus + phi_minus - 2*phi)*dx2_1;
//         }
//     }
//     nabla2 /= phi;

//     std::cout << "Analytical: " << nablaAnal << "   \t Numerical: " << nabla2 << "   \t Rel Diff: " << abs((nabla2-nablaAnal)/nabla2) << std::endl;
// }


std::vector<double> SimpleGaussian3D::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************
    double alpha = m_parameters[0], beta = m_parameters[1];
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    force[0] *= -2*alpha;
    force[1] *= -2*alpha;
    force[2] *= -2*alpha*beta;
    return force;
}

std::vector<double> SimpleGaussian3D::quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************
    double alpha = m_parameters[0], beta = m_parameters[1];
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    force[0] += step[0];
    force[0] *= -2*alpha;
    force[1] += step[1];
    force[1] *= -2*alpha;
    force[2] += step[2];
    force[2] *= -2*alpha*beta;
    return force;
}

double SimpleGaussian3D::phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    // Calculate (phi(new)/phi(old))**2
    auto pos = particles[index]->getPosition();
    double dr2=0, beta = m_parameters[1];
    dr2 += (2*pos[0]+step[0])*step[0]; 
    dr2 += (2*pos[1]+step[1])*step[1]; 
    dr2 += (2*pos[2]+step[2])*step[2]*beta; 
        // (pos[i]+step[i])*(pos[i]+step[i])-pos[i]*pos[i]
        // dr2 += (pos[i]+step[i])*(pos[i]+step[i])-pos[i]*pos[i];
    return exp(-2*m_parameters[0]*dr2);
}

std::vector<double> SimpleGaussian3D::getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles){
    double r2 = 0, z2 = 0, alpha = m_parameters[0], beta = m_parameters[1];
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        r2 += position[0]*position[0];
        r2 += position[1]*position[1];
        z2 += position[2]*position[2];
    }

    return std::vector<double>{-(r2+beta*z2), -alpha*z2};
}