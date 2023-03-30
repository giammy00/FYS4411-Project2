#pragma once
#include "sampler.h"
#include<fstream>
class SamplerFineTune: public Sampler {
public:
    SamplerFineTune(
    unsigned int numberOfParticles,
    unsigned int numberOfDimensions,
    int numberOfWFParams
    );

    std::ofstream m_outBinaryFile;
    void sample(bool acceptedStep, System* system) override ;
};
