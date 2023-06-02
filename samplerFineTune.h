#ifndef __SAMPLER_FINE_TUNE__
#define __SAMPLER_FINE_TUNE__
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
    SamplerFineTune(std::vector<std::unique_ptr< class SamplerFineTune  >  >  & samplers );
    void writeHistogram();
    std::ofstream m_outBinaryFile;
    std::ofstream m_outFileDistance;
    void sample(bool acceptedStep, System* system) override ;
    unsigned int *** m_position_histogram;
    unsigned int *** m_position_histogram2;
    unsigned int **** m_histograms;
    int m_nx=100;
    int m_ny=100;
    int m_nz=100;
    double  m_xMin, m_xMax, m_yMin , m_yMax, m_zMin,  m_zMax ; 
};
#endif