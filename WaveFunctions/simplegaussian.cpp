#include <memory>
#include <math.h>
#include <assert.h>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha > 0); // If alpha == 0 then the wavefunction doesn't go to zero
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    double phi = 1;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        for (unsigned int j = 0; j<position.size(); j++)
            phi *= exp(-1*position[j]*position[j]*m_parameters[0]);
    }
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
    /* The second derivative of exp(-alpha x*x) is exp(-alpha x*x)*(4*alpha*alpha*x*x - 2*alpha)
    */
    double sum = 0;
    for (unsigned int i = 0; i < particles.size(); i++){
        std::vector<double> position = particles[i]->getPosition();
        for (unsigned int j = 0; j < particles[i]->getNumberOfDimensions(); j++)
            sum += exp(-1*position[j]*position[j]*m_parameters[0]) * (4*position[j]*position[j]*m_parameters[0]*m_parameters[0] - 2*m_parameters[0]);
    }
    return sum;
}
