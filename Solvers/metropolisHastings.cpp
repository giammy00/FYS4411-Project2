#include <memory>
#include <vector>

#include "metropolisHastings.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"

#include <iostream>
using std::cout;
using std::endl;

MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool MetropolisHastings::step(
        double dt,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles)
{
    double temp;
    double sqrtdt = sqrt(dt);
    int index = m_rng->nextInt(0,particles.size()-1);
    auto step = waveFunction.quantumForce(particles, index);
    auto noise = std::vector<double>();
    for(unsigned int i=0; i<step.size(); i++){
        noise.push_back(m_rng->nextGaussian(0,sqrtdt));
        step[i] = step[i]*dt*0.5 + noise[i];
    }
    
    double hastingsArticle = waveFunction.partialHastingsArticle(particles, index, step);

    // argument of the greensfunction exponent: -(y-x-D*dt*F(x))^2 / 4*D*dt
    // difference: (-(y-x-D*dt*F(x))^2+(x-y-D*dt*F(y))^2) / 4*D*dt
    // difference: (-(step-D*dt*F(x))^2+(-step-D*dt*F(y))^2) / 4*D*dt
    // difference: (-(noise)^2+(-step-D*dt*F(y))^2) / 4*D*dt
    // difference: (-(noise)^2+(step+D*dt*F(y))^2) / 4*D*dt
    
    double greensDiff = 0;
    for (unsigned int i=0; i<noise.size(); i++)
        greensDiff -= noise[i]*noise[i];
    auto ForceY = waveFunction.quantumForceMoved(particles, index, step);
    for(unsigned int i=0; i<step.size(); i++){
        temp = ForceY[i]*dt*0.5 + step[i];
        greensDiff += temp*temp;
    }

    hastingsArticle *= exp(greensDiff/(2*dt));

    if(m_rng->nextDouble()<(hastingsArticle)){
        particles[index]->adjustPosition(step);
        return true;
    }
    return false;
}
