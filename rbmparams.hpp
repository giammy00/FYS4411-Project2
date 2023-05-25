#ifndef __RBM_PARAMS__
#define __RBM_PARAMS__
#include<vector>
#include"Math/random.h"
#include<string>

/*this class stores the parameters of a RBM in a std::vector m_allParams, and it implements
access to the parameters using conventional split into a, b , W, with proper indexing. 
m_a, m_b, m_W are pointers to the proper parts of m_allParams
In this way
(1) The code can interact with nlopt library, which requires gradient and parameters
to be stored in a std::vector 
(2) We can access the parameters easily using indexing, like m_a[i], and so on
*/
class RBMParams {
    public:
        RBMParams(unsigned int Nvisible,unsigned int Nhidden,  Random * rng);
        RBMParams(int Nvisible, int Nhidden, std::string filename);
        void setParams(std::vector<double> input);
        void saveParams(std::string filename);
        void readParams(std::string filename);
        double * m_a;
        double * m_b;
        double ** m_W;
        unsigned int m_Nvisible;
        unsigned int m_Nhidden;
        std::vector<double> m_allParams ;
};
#endif