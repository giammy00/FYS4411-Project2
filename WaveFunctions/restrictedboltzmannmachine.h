#pragma once
#include<iostream>
#include"wavefunction.h"
#include<memory>

class RestrictedBoltzmannMachine : public WaveFunction {

    public:
        RestrictedBoltzmannMachine(double sigma, int nVisible, int Nhidden, Random * rng);
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);

    private:
    //parameters of restricted boltzmann machine
        double m_sigma ;
        std::vector<double> m_a;
        std::vector<double> m_b;
        std::vector<std::vector<double>> m_W ;

};