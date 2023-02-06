#include "particle.h"
#include <assert.h>

Particle::Particle(const std::vector<double>& position) {
    m_numberOfDimensions = position.size();
    m_position = position;
}

void Particle::adjustPosition(double change, unsigned int dimension) {
    m_position.at(dimension) += change;
}
