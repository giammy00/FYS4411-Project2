#include<memory>
#include <cassert>
#include <iostream>
#include<cmath>
#include "interactingharmonicoscillator.hpp"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

InteractingHarmonicOscillator::InteractingHarmonicOscillator(double omega)
{
    assert(omega > 0);
    m_omega  = omega;
}

double InteractingHarmonicOscillator::computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles)
{

    /*Harmonic oscillator potential + coulomb term + kinetic energy.*/

    double rInverse, r2;
    unsigned int N = particles.size();
    // m = omega = 1
    double kineticEnergy   = waveFunction.computeDoubleDerivative(particles)*-0.5;
    double coulombEnergy=0;
    double potentialEnergy=0;
    //loop through particles

    //this can in principle be optimize, because we compute 1/r_{ij} more times than necessary.
    
    for (unsigned int i = 0; i < particles.size(); i++){
        auto pos = particles[i]->getPosition();
        //compute distance from all particles with j<i
        for (unsigned int j = 0; j<i; j++){
            auto pos2 = particles[j]->getPosition();
            r2 = 0;
            //pythagoras
            for (unsigned int k = 0; k<2; k++){
                r2 += (pos2[k] - pos[k])*(pos2[k] - pos[k]);
            }
            rInverse = 1.0/sqrt(r2);
            coulombEnergy+=rInverse;
        }
        potentialEnergy += pos[0] * pos[0];
        potentialEnergy += pos[1] * pos[1];

    }
    potentialEnergy*=0.5;
    //return energy per particle!
    return (kineticEnergy + potentialEnergy + coulombEnergy)/N;
}

