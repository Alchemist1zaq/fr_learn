#include "intFace.hpp"

#include "ele.hpp"
#include "flux.hpp"

void intFace::setupRightState(void) {
  // This is kinda messy, but avoids separate initialize function
  faceID_R = myInfo.IDR;
  relRot = myInfo.relRot;

  if (nDims == 2)
    nFptsR = eR->order + 1;
  else
    nFptsR = (eR->order + 1) * (eR->order + 1);

  /* --- Will have to introduce 'mortar' elements in the future [for
   * p-adaptation], but for now just force all faces to have same # of flux
   * points [order] --- */

  if (nFptsL != nFptsR)
    fatalError(
        "Mortar elements not yet implemented - must have nFptsL==nFptsR");

  FnR.resize(nFptsR);
  normR.setup(nFptsR, nDims);
  dAR.resize(nFptsR);

  /* --- Setup the L/R flux-point matching --- */
  fptR.resize(nFptsL);
  if (nDims == 2) {
    // For 1D faces [line segments] only - find first/last ID of fpts;
    // right element's are simply reversed
    fptStartR = (faceID_R * (nFptsR)) + nFptsR;
    fptEndR = (faceID_R * (nFptsR));

    int fpt = 0;
    for (int i = fptStartR - 1; i >= fptEndR; i--) {
      fptR[fpt] = i;
      fpt++;
    }
  } else if (nDims == 3) {
    // Only for quad tensor-product faces
    int order = sqrt(nFptsL) - 1;

    fptStartR = nFptsR * faceID_R;
    fptEndR = nFptsR * (faceID_R + 1);

    for (int i = 0; i < nFptsL; i++) {
      int ifpt = i % (order + 1);
      int jfpt = floor(i / (order + 1));
      switch (relRot) {
      case 0:
        fptR[i] = fptStartR + ifpt * (order + 1) + jfpt;
        break;
      case 1:
        fptR[i] = fptStartR + order - ifpt + jfpt * (order + 1);
        break;
      case 2:
        fptR[i] = fptEndR - 1 - (ifpt * (order + 1) + jfpt);
        break;
      case 3:
        fptR[i] = fptEndR - (order + 1) * (jfpt + 1) + ifpt;
        break;
      }
    }
  }

  if (params->viscous) {
    dUcR.setup(nFptsR, nFields);
  }
}

void intFace::getPointersRight(void) {
  // Get access to normal flux storage at right element [use look-up table to
  // get right fpt]
  for (int i = 0; i < nFptsL; i++) {
    FnR[i] = &(eR->Fn_fpts(fptR[i], 0));

    if (params->viscous)
      for (int k = 0; k < nFields; k++)
        dUcR(i, k) = &(eR->dUc_fpts(fptR[i], k));
  }
}

void intFace::getRightState(void) {
  // Get data from right element [order reversed to match left ele]
  for (int fpt = 0; fpt < nFptsL; fpt++) {
    for (int j = 0; j < nFields; j++) {
      UR(fpt, j) = (eR->U_fpts(fptR[fpt], j));
    }

    /* For dynamic grids (besides rigid translation), need to update
     * geometry-related data on every iteration, not just during setup */
    if (isNew_R || (params->motion != 0 && params->motion != 4)) {
      for (int dim = 0; dim < nDims; dim++) {
        normR(fpt, dim) = (eR->norm_fpts(fptR[fpt], dim));
      }
      dAR[fpt] = (eR->dA_fpts(fptR[fpt]));
    }
  }

  isNew_R = false;
}

void intFace::getRightGradient(void) {
  // Get data from right element [order reversed to match left ele]
  if (params->viscous) {
    for (int fpt = 0; fpt < nFptsL; fpt++) {
      for (int dim = 0; dim < nDims; dim++)
        for (int j = 0; j < nFields; j++)
          gradUR[fpt](dim, j) = (eR->dU_fpts(dim, fptR[fpt], j));
    }
  }
}

void intFace::setRightStateFlux(void) {
  for (int i = 0; i < nFptsR; i++)
    for (int j = 0; j < nFields; j++)
      FnR[i][j] = -Fn(i, j) * dAR[i]; // opposite normal direction
}

void intFace::setRightStateSolution(void) {
  for (int i = 0; i < nFptsR; i++)
    for (int j = 0; j < nFields; j++)
      *dUcR(i, j) = UC(i, j) - UR(i, j);
}

vector<double> intFace::computeWallForce() {
  // Not a wall boundary - return 0
  vector<double> force = {0, 0, 0, 0, 0, 0};
  return force;
}

vector<double> intFace::computeMassFlux() {
  // Not an inlet/outlet boundary - return 0
  vector<double> force(nFields);
  return force;
}
