#include <memory>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string>

#include "interactinggaussian3d.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

void testDoubleDerivative3d(
    std::vector<std::unique_ptr<class Particle>>& particles,
    class WaveFunction& waveFunction,
    double nabla2);
void slowEvaluate3d(std::vector<std::unique_ptr<class Particle>>& particles,
    std::vector<std::vector<double>>& m_distances,
    std::vector<double>& m_parameters, double phi_);

InteractingGaussian3D::InteractingGaussian3D(double alpha, double beta, double a)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_parameters.push_back(a);
}

void InteractingGaussian3D::InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles){
    double r2, dist, u_p;
    assert(particles[0]->getNumberOfDimensions() == 3);
    m_distances = std::vector<std::vector<double>>();
    for (unsigned int i = 0; i < particles.size(); i++){
        auto pos = particles[i]->getPosition();
        auto temp = std::vector<double>();
        for (unsigned int j = 0; j<i; j++){
            auto pos2 = particles[j]->getPosition();
            r2 = 0;
            for (unsigned int k = 0; k<3; k++){
                r2 += (pos2[k] - pos[k])*(pos2[k] - pos[k]);
            }
            temp.push_back(sqrt(r2));
        }
        r2 = 0;
        for (unsigned int k = 0; k<3; k++){
            r2 += pos[k] * pos[k];
        }
        temp.push_back(r2);
        m_distances.push_back(temp);
    }
    m_interForces = std::vector<double>(3*particles.size(),0);
    for (unsigned int k = 0; k < particles.size(); k++){
        auto pos = particles[k]->getPosition();
        for (unsigned int i = 0; i<particles.size(); i++){
            if (i == k) continue;
            auto pos2 = particles[i]->getPosition();
            dist = i<k? m_distances[k][i] : m_distances[i][k];
            u_p = uPrime_r(dist);
            for (unsigned int j = 0; j<3; j++){
                m_interForces[3*k+j] += (pos[j]-pos2[j]) * u_p;
            }
        }
    }
}

void InteractingGaussian3D::adjustPosition(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double> step){
    
    double r2, u_p;
    auto pos_old = particles[index]->getPosition();
    auto pos = particles[index]->getPosition();
    for (unsigned int j = 0; j<3; j++){
        m_interForces[3*index+j] = 0;
        pos[j] += step[j];
    }

    for (int j = 0; j<index; j++){ // Row
        auto pos2 = particles[j]->getPosition();
        u_p = uPrime_r(m_distances[index][j]);
        for (unsigned int k = 0; k<3; k++)
            m_interForces[3*j+k] -= (pos2[k] - pos_old[k])*u_p;
        r2 = 0;
        for (unsigned int k = 0; k<3; k++){
            r2 += (pos2[k] - pos[k])*(pos2[k] - pos[k]);
        }
        m_distances[index][j]=sqrt(r2);
        u_p = uPrime_r(m_distances[index][j]);
        for (unsigned int k = 0; k<3; k++){
            m_interForces[3*j+k] += (pos2[k] - pos[k])*u_p;
            m_interForces[3*index+k] += (pos[k]-pos2[k]) * u_p;
        }
    }
    r2 = 0;
    for (unsigned int k = 0; k<3; k++){
        r2 += pos[k] * pos[k];
    }
    m_distances[index][index] = r2;
    for (unsigned int j = index+1; j<particles.size(); j++){ // Column
        auto pos2 = particles[j]->getPosition();
        u_p = uPrime_r(m_distances[j][index]);
        for (unsigned int k = 0; k<3; k++)
            m_interForces[3*j+k] -= (pos2[k] - pos_old[k])*u_p;
        r2 = 0;
        for (unsigned int k = 0; k<3; k++){
            r2 += (pos2[k] - pos[k])*(pos2[k] - pos[k]);
        }
        m_distances[j][index]=sqrt(r2);
        u_p = uPrime_r(m_distances[j][index]);
        for (unsigned int k = 0; k<3; k++){
            m_interForces[3*j+k] += (pos2[k] - pos[k])*u_p;
            m_interForces[3*index+k] += (pos[k]-pos2[k]) * u_p;
        }
    }
}

double InteractingGaussian3D::uPrime_r(double r){
    // *************** [Derivative of log(1.-(m_parameters[2]/r))]/r ***************
    double a = m_parameters[2];
    return a/(r*r*abs(a-r));
}

double InteractingGaussian3D::uDoublePrime(double r){
    // *************** Second derivative of log(1.-(m_parameters[2]/r)) ***************
    double a = m_parameters[2];
    return (a*(a-2*r))/(r*r*(a-r)*(a-r));
    // return 1/(r*r)-1/((m_parameters[2]-r)*(m_parameters[2]-r));
}

double InteractingGaussian3D::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    // Returns Phi, not Phi^2
    double r2 = 0, a = m_parameters[2];
    for (unsigned int i = 0; i < particles.size(); i++){
        r2 += m_distances[i][i];
    }
    double phi = exp(-1*r2*m_parameters[0]);
    for (unsigned int i = 0; i < particles.size(); i++){
        for (unsigned int j = i+1; j<particles.size(); j++){
            phi *= std::max(1-a/m_distances[j][i],0.);
        }
    }

    // TESTING:
    // slowEvaluate(particles, m_distances, m_parameters, phi);

    return phi;
}

void slowEvaluate3d(std::vector<std::unique_ptr<class Particle>>& particles, std::vector<std::vector<double>>& m_distances, std::vector<double>& m_parameters, double phi_) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    // Returns Phi, not Phi^2
    double r2 = 0;
    double diffSum=0;
    int count=0;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> pos = particles[i]->getPosition();
        for (unsigned int j = 0; j<3; j++)
            r2 += pos[j]*pos[j];
    }
    double phi = exp(-1*r2*m_parameters[0]);
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> pos = particles[i]->getPosition();
        for (unsigned int j = i+1; j<particles.size(); j++){
            std::vector<double> pos2 = particles[j]->getPosition();
            r2 = 0;
            for (unsigned int k = 0; k<3; k++)
                r2 += (pos[k]-pos2[k])*(pos[k]-pos2[k]);
            phi *= std::max(1-m_parameters[2]/sqrt(r2),0.);
            count++;
            diffSum += abs((sqrt(r2)-m_distances[j][i])/m_distances[j][i]);
            // std::cout << sqrt(r2)-m_distances[j][i] << "\t";
        }
        // std::cout << std::endl;
    }
    std::cout << "Rel phi diff: " << abs((phi-phi_)/phi) << "\t Avg relative r diff: " << diffSum/count << std::endl;
}

double InteractingGaussian3D::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) {
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
    double r2 = 0, alpha = m_parameters[0];
    for (unsigned int i = 0; i < particles.size(); i++){
        r2 += m_distances[i][i];
    }
    int n = particles.size();
    double nabla2 = 4*alpha*alpha*r2 - 6*n*alpha;

    // double sum over all particles
    double dist, u_p;
    std::vector<double> repulsion, position;
    for (unsigned int k = 0; k < particles.size(); k++){
        for (unsigned int i = 0; i < k; i++){           
            dist = m_distances[k][i];
            u_p = uPrime_r(dist);
            nabla2 += 2*uDoublePrime(dist) + 4*u_p;
        }
        position = particles[k]->getPosition();
        for (unsigned int j = 0; j < 3; j++){
            // nabla phi = -2*alpha*r
            nabla2 += -4 * alpha * position[j] * m_interForces[3*k+j];
            nabla2 += m_interForces[3*k+j] * m_interForces[3*k+j];
        }
    }

    // TESTING:
    // testDoubleDerivative(particles, *this, nabla2);
    // slowDoubleDerivative(particles, nabla2);

    return nabla2;
}

void InteractingGaussian3D::slowDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles, double nabla2_) {
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
        for (unsigned int j = 0; j < 3; j++)
            r2 += position[j]*position[j];
    }
    int n = particles.size();
    double nabla2 = 4*m_parameters[0]*m_parameters[0]*r2 - 3*n*m_parameters[0];
    // int n = particles.size() * 3;
    // double nabla2 = 4*m_parameters[0]*m_parameters[0]*r2 - 2*n*m_parameters[0];

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
            for (unsigned int j = 0; j < 3; j++){
                diff = position[j] - position2[j];
                r2 += diff*diff;
            }
            dist = sqrt(r2);
            u_p = uPrime_r(dist);
            for (unsigned int j = 0; j < 3; j++){
                diff = position[j] - position2[j];
                repulsion[j] += diff*u_p;
            }
            nabla2 += uDoublePrime(dist) + 2*u_p;
        }
        for (unsigned int j = 0; j < 3; j++){
            // nabla phi = -2*alpha*r
            nabla2 += -4 * m_parameters[0] * position[j] * repulsion[j];
            nabla2 += repulsion[j] * repulsion[j];
        }
    }

    std::cout << "Nabla^2 diff: " << abs((nabla2-nabla2_)/nabla2) << std::endl;
}
    
void testDoubleDerivative3d(
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
        for (unsigned int j = 0; j < 3; j++){
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


std::vector<double> InteractingGaussian3D::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    double alpha = m_parameters[0];
    for (unsigned int i=0; i<3; i++){
        force[i] *= -2*alpha;
        force[i] += m_interForces[3*index+i];
    }
    // TESTING:
    // slowQuantumForce(particles, index, force);
    return force;
}

void InteractingGaussian3D::slowQuantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& force_){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************
    double r2, temp;
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    for (unsigned int i=0; i<3; i++){
        force[i] *= -2*m_parameters[0];
    }
    for (unsigned int i = 0; i<particles.size(); i++){
        if ((int)i == index) continue;
        auto pos2 = particles[i]->getPosition();
        auto relPos = std::vector<double>();
        r2 = 0;
        for (unsigned int k = 0; k<3; k++){
            relPos.push_back(pos[k] - pos2[k]);
            r2 += relPos[k]*relPos[k];
        }
        temp = sqrt(r2);
        temp = uPrime_r(temp);
        for (unsigned int k = 0; k<3; k++){
            force[k] += relPos[k] * temp;
        }
    }
    double diff = 0;
    for (unsigned int k = 0; k<3; k++)
        diff += (force[k]-force_[k])*(force[k]-force_[k]);
    std::cout << "force diff: " << diff << std::endl;
}

std::vector<double> InteractingGaussian3D::quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    //***************WE RETURN d/dx(phi)/phi NOT d/dx(phi)*********************

    double r2, temp, alpha = m_parameters[0];
    auto pos = particles[index]->getPosition();
    auto force = std::vector<double>(pos);
    for (unsigned int i=0; i<3; i++){
        force[i] += step[i];
        force[i] *= -2*alpha;
    }
    for (unsigned int i = 0; i<particles.size(); i++){
        if ((int)i == index) continue;
        auto pos2 = particles[i]->getPosition();
        auto relPos = std::vector<double>{0,0,0};
        r2 = 0;
        for (unsigned int k = 0; k<3; k++){
            relPos[k] = (pos[k] + step[k] - pos2[k]);
            r2 += relPos[k]*relPos[k];
        }
        temp = sqrt(r2);
        temp = uPrime_r(temp);
        for (unsigned int k = 0; k<3; k++){
            force[k] += relPos[k] * temp;
        }
    }
    return force;
}

double InteractingGaussian3D::phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step){
    // Calculate (phi(new)/phi(old))**2
    double r2 = 0, r2new = 0, diff, Jik, JikNew, a = m_parameters[2];
    auto pos = particles[index]->getPosition();
    for (unsigned int i = 0; i<3; i++)
        r2 += (pos[i]+step[i])*(pos[i]+step[i])-pos[i]*pos[i];
    double phi = exp(-2*r2*m_parameters[0]);
    for (unsigned int i = 0; i<particles.size(); i++){
        if ((int)i == index) continue;
        auto pos2 = particles[i]->getPosition();
        r2 = 0;
        r2new = 0;
        for (unsigned int k = 0; k<3; k++){
            diff = pos[k] - pos2[k];
            r2 += diff*diff;
            r2new += (diff+step[k])*(diff+step[k]);
        }
        Jik = 1-a/sqrt(r2);
        if (Jik<0)
            phi *= 1E6; // arbitrary constant to motivate particles that overlap to move
        else{
            JikNew = std::max(1-a/sqrt(r2new),0.);
            phi *= JikNew*JikNew/(Jik*Jik);
        }
    }
    return phi;
}