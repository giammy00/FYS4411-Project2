#include <memory>
#include <iostream>
#include <cassert>

#include "initialstate.h"
#include "Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng,
            double lengthScale
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    for (unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (unsigned int j=0; j < numberOfDimensions; j++) {
            position.push_back(rng.nextDouble(-lengthScale, lengthScale));
        }
        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}

std::vector<std::unique_ptr<Particle>> setupRandomGaussianInitialState(
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng,
            double lengthScale
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    for (unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (unsigned int j=0; j < numberOfDimensions; j++) {
            position.push_back(rng.nextGaussian(0,lengthScale));
        }
        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}
