#pragma once
#include "sampler.h"
#include<fstream>
class SamplerFineTune: public Sampler {
public:
//this sampler is a simplified version of Sampler which does not sample the quantities used to compute the gradient. 
    SamplerFineTune(
    unsigned int numberOfParticles,
    unsigned int numberOfDimensions,
    int numberOfWFParams
    );

    std::ofstream m_outBinaryFile;
    void sample(bool acceptedStep, System* system) override ;
    unsigned int *** m_position_histogram;
};
