#pragma once
#include "sampler.h"
#include<fstream>
class SamplerFineTune: public Sampler {
public:
    SamplerFineTune(
    unsigned int numberOfParticles,
    unsigned int numberOfDimensions,
    int numberOfWFParams,
    int thread_number
    );

    std::ofstream m_outBinaryFile;
    void SamplerFineTune::sample(bool acceptedStep, System* system) ;
};
