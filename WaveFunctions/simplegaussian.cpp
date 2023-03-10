#include <memory>
#include <math.h>
#include <assert.h>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

SimpleGaussian::SimpleGaussian(double alpha, std::vector<std::unique_ptr<class Particle>>& particles)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    r2 = 0;
    //store sum of positions squared, useful in many computations
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        for (unsigned int j = 0; j<position.size(); j++)
            r2 += position[j]*position[j];
    }
    
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    // Returns Phi, not Phi^2
    // double r2 = 0;
    // for (unsigned int i = 0; i < particles.size(); i++){
    //     std::vector<double> position = particles[i]->getPosition();
    //     for (unsigned int j = 0; j<position.size(); j++)
    //         r2 += position[j]*position[j];
    // }
    double phi = exp(-1*r2*m_parameters[0]);
    return phi;
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
    // double r2 = 0;
    // for (unsigned int i = 0; i < particles.size(); i++){
    //     std::vector<double> position = particles[i]->getPosition();
    //     for (unsigned int j = 0; j < particles[i]->getNumberOfDimensions(); j++)
    //         r2 += position[j]*position[j];
    // }
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

std::vector<double> SimpleGaussian::getGradientTerms( double Elocal ){
    // get the term necessary to the computation of cumulative gradient i.e. (1/psi)(d/d\alpha)psi = (d/dalpha)log(psi)
    double dLogPhi_dAlpha = -r2;
    double dLogPhiTimesEL = dLogPhi_dAlpha*Elocal;
    std::vector<double> gradTerms = std::vector<double>{dLogPhi_dAlpha, dLogPhiTimesEL};
    return gradTerms;
    }

void SimpleGaussian::updateCachedVariables(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<double>& step){
    //this function helps updating some stored quantities which come up again and again in various computations, to spare computational costs.
    //in the case of harmonic oscillator the only quantity is: r2=sum_i^Nparticles ( r_i^2 )
    //which changes by 2*step*(1+step), when we move EXACTLY one particle by a quantity step (regardless of which particle it is)
    for(unsigned int i=0; i<step.size(); i++){
        r2+=2*step[i]*(1+step[i]);
    }

    
 }