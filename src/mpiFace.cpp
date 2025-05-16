

#include "../include/mpiFace.hpp"

#include "../include/ele.hpp"
#include "../include/flux.hpp"

void mpiFace::setupRightState(void) {}

void mpiFace::finishRightSetup(void) {}

void mpiFace::getPointersRight() {
  // Do nothing
}

void mpiFace::communicate(void) { getLeftState(); }

void mpiFace::communicateGrad(void) {}

void mpiFace::getRightState(void) {}

void mpiFace::getRightGradient(void) {}

void mpiFace::setRightStateFlux(void) {}

void mpiFace::setRightStateSolution(void) {}

vector<double> mpiFace::computeWallForce() {
  // Not a wall boundary - return 0
  vector<double> force = {0, 0, 0, 0, 0, 0};
  return force;
}

vector<double> mpiFace::computeMassFlux() {
  // Not an inlet/outlet boundary - return 0
  vector<double> force(nFields);
  return force;
}
