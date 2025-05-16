
#pragma once

class face;
class overComm;

#include "face.hpp"

#include "overComm.hpp"

class overFace : public face {
public:
  //! Overset communicator object where the face's data will be found
  overComm *OComm;

  //! Setup arrays to handle getting/setting data from Solver
  void setupRightState(void);

  /*! Get pointer access to right ele's data */
  void getPointersRight(void);

  //! Get the interpolated overset solution data from the Solver
  void getRightState(void);

  //! Get the interpolated overset gradient data from the Solver
  void getRightGradient(void);

  //! Do nothing [right state is non-existant]
  void setRightStateFlux(void);

  //! Do nothing [right state is non-existant]
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

  //! Do nothing [not an inlet/outlet boundary]
  vector<double> computeMassFlux(void);

  //! Return the physical position of the face's flux points
  vector<point> getPosFpts(void);

  //! Return the outward unit normal at the face's flux points
  vector<point> getNormFpts();

  //! Override normal version when using flux-interp method
  void rusanovFlux(void);

  int fptOffset; //! Offset within Solver's mesh block-global interp point list
  vector<point> posFpts; //! Physical locations of left ele's flux points

private:
};
