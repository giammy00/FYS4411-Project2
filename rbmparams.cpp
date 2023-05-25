#include"rbmparams.hpp"
#include<fstream>
#include<iostream>
RBMParams::RBMParams(unsigned int Nvisible,unsigned int Nhidden,  Random * rng){
    //random initializer for trainable parameters of the restricted boltzmann machine 
    
    m_allParams.reserve(m_Nvisible+m_Nhidden+m_Nvisible*m_Nhidden);
    m_Nvisible =Nvisible;
    m_Nhidden = Nhidden;
    double tmp;
    //these factors are used to scale the initial random parameters. Arbitrary choice.
    double norm_a = 1.0/Nvisible;
    double norm_b = 1.0/Nhidden;
    double norm_w = norm_a*norm_b;
    //init weights a
    for(size_t i = 0; i<Nvisible; i++){
        tmp = rng->nextGaussian(0, norm_a);
        m_allParams.push_back(tmp);
    }
    m_a = m_allParams.data();

    //init weights b
    for(size_t i = 0; i<Nhidden; i++){
        tmp = rng->nextGaussian(0, norm_b);
        m_allParams.push_back(tmp);
        }
    //init pointers. 
    m_b = m_a+m_Nvisible;

    m_W = (double **) malloc(Nvisible*sizeof(double*));
    //init weights W
    for(size_t i = 0; i<Nvisible; i++){
        std::vector<double> tmp_arr(Nhidden);
        m_W[i] = m_b+Nhidden+i*Nhidden;
        //fill tmp array with random numbers
        for(size_t j =0; j<Nhidden; j++){
            tmp = rng->nextGaussian(0, norm_w);
            tmp_arr[j] = tmp;
            m_allParams.push_back(tmp);
        }
    }
    return ; 
}
RBMParams::RBMParams(int Nvisible, int Nhidden, std::string filename){
    m_Nvisible = Nvisible;
    m_Nhidden = Nhidden;
    readParams(filename);
    return ; 
}

void RBMParams::setParams( std::vector<double> params ){
    //sets the parameters given a flattened vector
    //nb order matters:
    //[visible_params, hidden_params, WeightMatrix(flattened)]

    //copy weights a
    unsigned int offset = 0;
    for(size_t i = 0; i<m_Nvisible; i++){
        m_a[i] = params[offset+i] ;
        m_allParams[offset+i] = params[offset+i];
    }
    offset+=m_Nvisible;
    //copy weights b
    for(size_t i = 0; i<m_Nhidden; i++){
        m_b[i]= params[offset+i];   
        m_allParams[offset+i] = params[offset+i];
        }
    //copy weights W
    offset+=m_Nhidden;
    for(size_t i = 0; i<m_Nvisible; i++){
        //fill tmp array with random numbers
        for(size_t j =0; j<m_Nhidden; j++){
            m_W[i][j]= params[offset+i*m_Nhidden+j];
            m_allParams[offset+i] = params[offset+i];
        }
    }
    return ; 
}

void RBMParams::saveParams(std::string filename){
    //save the parameters to a binary file to be read later
    std::ofstream outBinaryFile(filename, std::ios::binary);
    outBinaryFile.write(reinterpret_cast<const char*>(m_allParams.data()), m_allParams.size() * sizeof(double));
    outBinaryFile.close();
    return ; 
}

void RBMParams::readParams(std::string filename){
    //read the parameters saved in filename, into the std::vector m_allParams
    //also set the pointers
    std::ifstream file(filename, std::ios::binary);
    if (file.is_open()) {
        // Read the elements into the vector
        unsigned int size = m_Nvisible+m_Nhidden+m_Nvisible*m_Nhidden;
        m_allParams.reserve(size);
        file.read(reinterpret_cast<char*>(m_allParams.data()), size*sizeof(double) );
        file.close();
        m_a = m_allParams.data();
        m_b = m_a + m_Nvisible;
        m_W = (double **) malloc(m_Nvisible*sizeof(double*))    ;
        for (size_t i=0; i<m_Nvisible; i++){
            m_W[i] = m_b+m_Nhidden+i*m_Nhidden;
        }    
        std::cout << "Successfully read parameters." << std::endl;
    } else {
        std::cout << "Failed to open " << filename << ". Terminating." << std::endl;
        exit(1);
    }
    return ; 
}