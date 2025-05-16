
#pragma once

#include <vector>

#include "face.hpp"
#include "matrix.hpp"

class face;

class mpiFace : public face {
public:
  //! Get/send information from/to the opposite processor
  void setupRightState(void);

  //! Wait to receive # of flux points from opposite proccessor
  void finishRightSetup(void);

  /*! Get pointer access to right ele's data */
  void getPointersRight(void);

  //! Receive the right-state data from the opposite processor
  void getRightState(void);

  //! For viscous cases, receive the solution gradient from the opposite
  //! processor
  void getRightGradient(void);

  //! Do nothing [handled sparately via comminicate()]
  void setRightStateFlux(void);

  //! Do nothing [handled sparately via comminicate()]
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

  //! Do nothing [not an inlet/outlet boundary]
  vector<double> computeMassFlux(void);

  //! Send the right-state data across the processor boundary using MPI
  void communicate(void);

  //! Send the right-state gradient data across the processor boundary using MPI
  void communicateGrad();

  int procL; //! Processor ID on left  [this face]
  int procR; //! Processor ID on right [opposite face]
  int IDR;   //! Local face ID of face on right processor

private:
  int faceID_R; //! Right element's element-local face ID
  int relRot;   //! Relative rotation of right element's face (for 3D)
  int fptStartR, fptEndR;
  vector<int> fptR; //! Indices of flux points on right face

  /* --- Storage for all solution/geometry data at flux points [right state] ---
   */
  // vector<matrix<double>> FR;   //! Flux array [nFpts, nDims, nFields]
  matrix<double> normR;   //! Unit outward normal at flux points  [nFpts, nDims]
  vector<double> dAR;     //! Local face-area equivalent at flux points
  vector<double> detJacR; //! Determinant of transformation Jacobian

  matrix<double> bufUR;       //! Incoming buffer for receving UR
  Array<double, 3> bufGradUR; //! Incoming buffer for receving gradUR
  Array<double, 3> bufGradUL; //! !! TEMP HACK !! Outgoing buffer for sending UL
};
