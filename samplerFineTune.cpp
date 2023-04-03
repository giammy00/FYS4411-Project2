#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <omp.h>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include "samplerFineTune.h"
#include "utils.h"
using std::cout;
using std::endl;


SamplerFineTune::SamplerFineTune( unsigned int numberOfParticles,
    unsigned int numberOfDimensions,
    int ) : Sampler::Sampler(numberOfParticles, numberOfDimensions, 0)
{
    int thread_number = omp_get_thread_num();
    std::string fname_thread = "./Outputs/sampledEnergies_" + std::to_string(numberOfParticles) + "_" + std::to_string(thread_number) + ".bin";
    m_outBinaryFile.open(fname_thread, std::ios::binary);
	
    int nx=100;
    int ny=100;
    int nz=100;
    m_position_histogram = init_3d_array(nx, ny, nz);
}
void SamplerFineTune::sample(bool acceptedStep, System* system) {
    /*sample all the interesting things 
     */
    double localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += localEnergy * localEnergy;

    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;

    //write sampled energy to a file
    m_outBinaryFile.write(reinterpret_cast<const char*>(&localEnergy), sizeof(double));

}
