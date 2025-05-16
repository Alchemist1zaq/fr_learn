
#pragma once

class face;

#include "face.hpp"

class intFace : public face {
public:
  //! Setup arrays to handle getting/setting data at right element
  void setupRightState(void);

  /*! Get pointer access to right ele's data */
  void getPointersRight(void);

  //! Get data from the right element
  void getRightState(void);

  //! For viscous cases, get the solution gradient from the right element
  void getRightGradient(void);

  //! Put the calculated interface flux into the right element's memory
  void setRightStateFlux(void);

  //! Put the common solution into the right element's memory (viscous cases)
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

  //! Do nothing [not an inlet/outlet boundary]
  vector<double> computeMassFlux(void);

private:
  int faceID_R; //! Right element's face ID
  int relRot;   //! Relative rotation of right element's face (for 3D)
  int fptStartR, fptEndR;
  vector<int> fptR; //! Indices of flux points in right element

  bool isNew_R = true; //! Flag for initialization (esp. due to unblanking)

  /* --- Storage for all solution/geometry data at flux points [right state] ---
   */
  vector<matrix<double>> FR; //! Flux array [nFpts, nDims, nFields]
  vector<double *> FnR; //! Common normal flux for right ele [in ele's memory]
  Array<double *, 2>
      dUcR; //! Common solution for left ele (in ele's memory)  [nFpts, nFields]
  matrix<double> normR; //! Unit outward normal at flux points  [nFpts, nDims]
  vector<double> dAR;   //! Local face-area equivalent at flux points
};
