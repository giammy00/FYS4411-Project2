#include <memory>
#include <math.h>
#include <assert.h>

#include "interactinggaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

InteractingGaussian::InteractingGaussian(double alpha, double beta, double a)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(a);
}

double InteractingGaussian::uPrime(double r){
    // *************** Derivative of log(1.-(m_parameters[2]/r)) ***************
    return m_parameters[2]/(r*abs(m_parameters[2]-r));
}

double InteractingGaussian::uDoublePrime(double r){
    // *************** Second derivative of log(1.-(m_parameters[2]/r)) ***************
    return (m_parameters[2]*(m_parameters[2]-2*r))/(r*r*(m_parameters[2]-r)*(m_parameters[2]-r));
    // return 1/(r*r)-1/((m_parameters[2]-r)*(m_parameters[2]-r));
}

double InteractingGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    // Returns Phi, not Phi^2
    double r2 = 0;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> pos = particles[i]->getPosition();
        for (unsigned int j = 0; j<pos.size(); j++)
            r2 += pos[j]*pos[j];
    }
    double phi = exp(-1*r2*m_parameters[0]);
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> pos = particles[i]->getPosition();
        for (unsigned int j = i+1; j<particles.size(); j++){
            std::vector<double> pos2 = particles[j]->getPosition();
            r2 = 0;
            for (unsigned int k = 0; k<pos.size(); k++)
                r2 += (pos[k]-pos2[k])*(pos[k]-pos2[k]);
            phi *= std::max(1-m_parameters[2]/sqrt(r2),0.);
        }
    }
    return phi;
}

double InteractingGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
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

    return nabla2;
    // This line always closes a multiline comment */
    
    
    
    /* Numerical calculation
    double nabla2 = 0, phi, phi_plus, phi_minus;
    // double phi, phi_plus, phi_minus;
    // nabla2 = 0;
    const double dx = 1e-5, dx2_1 = 1/(dx*dx); // dx2_1 = 1/(dx*dx)
    assert(abs(dx2_1 - 1/(dx*dx))<1);
    phi = evaluate(particles);
    for (unsigned int i = 0; i < particles.size(); i++){
        for (unsigned int j = 0; j < particles[i]->getNumberOfDimensions(); j++){
            particles[i]->adjustPosition(dx, j);
            phi_plus = evaluate(particles);
            particles[i]->adjustPosition(-2*dx, j);
            phi_minus = evaluate(particles);
            particles[i]->adjustPosition(dx, j);
            
            nabla2 += (phi_plus + phi_minus - 2*phi)*dx2_1;
        }
    }

    return nabla2/phi;
    // This line always closes a multiline comment */
}


std::vector<double> InteractingGaussian::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************
    double r2, temp;
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    for (unsigned int i=0; i<force.size(); i++){
        force[i] *= -2*m_parameters[0];
    }
    for (unsigned int i = 0; i<particles.size(); i++){
        if ((int)i == index) continue;
        auto pos2 = particles[i]->getPosition();
        auto relPos = std::vector<double>();
        r2 = 0;
        for (unsigned int k = 0; k<pos.size(); k++){
            relPos.push_back(pos[k] - pos2[k]);
            r2 += relPos[k]*relPos[k];
        }
        temp = sqrt(r2);
        temp = uPrime(temp)/temp;
        for (unsigned int k = 0; k<pos.size(); k++){
            force[k] += relPos[k] * temp;
        }
    }
    return force;
}

std::vector<double> InteractingGaussian::quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************

    double r2, temp;
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    for (unsigned int i=0; i<force.size(); i++){
        force[i] += step[i];
        force[i] *= -2*m_parameters[0];
    }
    for (unsigned int i = 0; i<particles.size(); i++){
        if ((int)i == index) continue;
        auto pos2 = particles[i]->getPosition();
        auto relPos = std::vector<double>();
        r2 = 0;
        for (unsigned int k = 0; k<pos.size(); k++){
            relPos.push_back(pos[k] + step[k] - pos2[k]);
            r2 += relPos[k]*relPos[k];
        }
        temp = sqrt(r2);
        temp = uPrime(temp)/temp;
        for (unsigned int k = 0; k<pos.size(); k++){
            force[k] += relPos[k] * temp;
        }
    }
    return force;
}

double InteractingGaussian::phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    // Calculate (phi(new)/phi(old))**2
    double r2 = 0, r2new = 0, diff, Jik, JikNew;
    std::vector<double> pos = particles[index]->getPosition();
    for (unsigned int i = 0; i<pos.size(); i++)
        r2 += (pos[i]+step[i])*(pos[i]+step[i])-pos[i]*pos[i];
    double phi = exp(-2*r2*m_parameters[0]);
    for (unsigned int i = 0; i<particles.size(); i++){
        if ((int)i == index) continue;
        std::vector<double> pos2 = particles[i]->getPosition();
        r2 = 0;
        r2new = 0;
        for (unsigned int k = 0; k<pos.size(); k++){
            diff = pos[k] - pos2[k];
            r2 += diff*diff;
            r2new += (diff+step[i])*(diff+step[k]);
        }
        Jik = 1-m_parameters[2]/sqrt(r2);
        if (Jik<0)
            phi *= 1E6; // arbitrary constant to motivate particles that overlap to move
        else{
            JikNew = std::max(1-m_parameters[2]/sqrt(r2new),0.);
            phi *= JikNew*JikNew/(Jik*Jik);
        }
    }
    return phi;
}