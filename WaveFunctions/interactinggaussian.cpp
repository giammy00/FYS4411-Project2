#include <memory>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string>

#include "interactinggaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

void testDoubleDerivative(
    std::vector<std::unique_ptr<class Particle>>& particles,
    class WaveFunction& waveFunction,
    double nabla2);

InteractingGaussian::InteractingGaussian(double alpha, double beta, double a)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(a);
}

void InteractingGaussian::InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles){
    double r2;
    distances = std::vector<std::vector<double>>();
    for (unsigned int i = 0; i < particles.size(); i++){
        auto pos = particles[i]->getPosition();
        auto temp = std::vector<double>();
        for (unsigned int j = 0; j<i; j++){
            auto pos2 = particles[j]->getPosition();
            for (unsigned int k = 0; k<pos.size(); k++){
                r2 += (pos2[k] - pos[k])*(pos2[k] - pos[k]);
            }
            temp.push_back(sqrt(r2));
        }
        distances.push_back(temp);
    }
}

double InteractingGaussian::uPrime_r(double r){
    // *************** [Derivative of log(1.-(m_parameters[2]/r))]/r ***************
    return m_parameters[2]/(r*r*abs(m_parameters[2]-r));
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
    
    
    // Non-interacting part
    double r2 = 0;
    for (unsigned int i = 0; i < particles.size(); i++){
        auto position = particles[i]->getPosition();
        for (unsigned int j = 0; j < particles[i]->getNumberOfDimensions(); j++)
            r2 += position[j]*position[j];
    }
    int n = particles.size() * particles[0]->getNumberOfDimensions();
    double nabla2 = 4*m_parameters[0]*m_parameters[0]*r2 - 2*n*m_parameters[0];

    // double sum over all particles
    double diff, dist, u_p;
    std::vector<double> repulsion, position;
    for (unsigned int k = 0; k < particles.size(); k++){
        repulsion = std::vector<double>{0.,0.,0.};
        position = particles[k]->getPosition();
        for (unsigned int i = 0; i < particles.size(); i++){
            if (i==k) continue;
            auto position2 = particles[i]->getPosition();
            r2 = 0;
            for (unsigned int j = 0; j < position.size(); j++){
                diff = position[j] - position2[j];
                r2 += diff*diff;
            }
            dist = sqrt(r2);
            u_p = uPrime_r(dist);
            for (unsigned int j = 0; j < position.size(); j++){
                diff = position[j] - position2[j];
                repulsion[j] += diff*u_p;
            }
            nabla2 += uDoublePrime(dist) + 2*u_p;
        }
        for (unsigned int j = 0; j < position.size(); j++){
            // nabla phi = -2*alpha*r
            nabla2 += -4 * m_parameters[0] * position[j] * repulsion[j];
            nabla2 += repulsion[j] * repulsion[j];
        }

    }

    // TESTING:
    // testDoubleDerivative(particles, *this, nabla2);

    return nabla2;
}
    
void testDoubleDerivative(
    std::vector<std::unique_ptr<class Particle>>& particles,
    class WaveFunction& waveFunction,
    double nablaAnal){
    
    //* Numerical calculation
    double nabla2 = 0, phi, phi_plus, phi_minus;
    // double phi, phi_plus, phi_minus;
    // nabla2 = 0;
    const double dx = 1e-5, dx2_1 = 1/(dx*dx);
    phi = waveFunction.evaluate(particles);
    for (unsigned int i = 0; i < particles.size(); i++){
        for (unsigned int j = 0; j < particles[i]->getNumberOfDimensions(); j++){
            particles[i]->adjustPosition(dx, j);
            phi_plus = waveFunction.evaluate(particles);
            particles[i]->adjustPosition(-2*dx, j);
            phi_minus = waveFunction.evaluate(particles);
            particles[i]->adjustPosition(dx, j);
            
            nabla2 += (phi_plus + phi_minus - 2*phi)*dx2_1;
        }
    }
    nabla2 /= phi;

    std::cout << "Analytical: " << nablaAnal << "   \t Numerical: " << nabla2 << "   \t Rel Diff: " << abs((nabla2-nablaAnal)/nabla2) << std::endl;
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
        temp = uPrime_r(temp);
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
        temp = uPrime_r(temp);
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