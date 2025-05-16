
#pragma once

#include "global.hpp"

#include <map>
#include <set>
#include <unordered_set>
#include <vector>

class oper;
class overFace;

#include "ele.hpp"
#include "input.hpp"
#include "operators.hpp"
#include "overFace.hpp"
#include "superMesh.hpp"

class overComm {
public:
  overComm();

  int nDims;
  int nGrids;
  int nprocPerGrid;
  int gridID;
  int gridRank;
  int rank;
  int nproc;

  vector<int> gridIdList; //! Grid ID for each MPI rank

  int nFields;

  input *params; //! Simulation input parameters

  /* --- Variables for Exchanging Data at Overset Faces --- */

  vector<int>
      nPts_rank; //! Number of fringe points for each rank of current grid
  vector<vector<int>> foundPts; //! IDs of receptor points from each grid which
                                //! were found to lie within current grid
  vector<vector<int>>
      foundRank; //! gridRank of this process for each found point (for benefit
                 //! of other processes; probably not needed)
  vector<vector<int>>
      foundEles; //! Ele ID which each matched point was found to lie within
  vector<vector<point>> foundLocs; //! Reference location within donor ele of
                                   //! each matched receptor point
  vector<vector<point>> foundNorm; //! For corrected-flux method: Outward unit
                                   //! normal at fringe boundary point

  vector<vector<matrix<double>>>
      foundJaco; //! For static cases with flux-interp, store transformation
                 //! matrix for each point
  vector<vector<double>> foundDetJac; //! For static cases with flux-interp,
                                      //! store Jacobian det. for each point

  int nOverPts; //! Number of overset (receptor) points on this grid
  matrix<double>
      overPts; //! Physical locations of the receptor points on this grid
  matrix<double> overNorm; //! Outward unit normals at fringe-boundary points on
                           //! this grid [corrected-flux interp method]
  vector<int>
      nPtsRecv; //! Number of points incoming from each grid (across interComm)
  vector<int>
      nPtsSend; //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>> recvPts; //! Point IDs which will be received from each
                               //! grid (across interComm) (counter to foundPts)

  matrix<double> U_in; //! Data received from other grid(s)
  vector<matrix<double>>
      U_out; //! Interpolated data being sent to other grid(s)

  matrix<double> gradU_in; //! Gradient data received from other grid(s)
  vector<matrix<double>>
      gradU_out; //! Interpolated gradient data being sent to other grid(s)

  /* --- Variables for Exchanging Data on Unblanked Cells --- */

  vector<int>
      nCells_rank; //! Number of unblanked cells for each rank of current grid
  vector<vector<int>> foundCells;      //! IDs of unblanked cells from each grid
                                       //! which were found to overlap this grid
  vector<matrix<int>> foundCellDonors; //! Donor cells on this grid for each
                                       //! found unblanked cell
  vector<vector<int>> foundCellNDonors; //! Number of donor cells on this grid
                                        //! for each found unblanked cell

  int nUnblanks;          //! Number of unblank cells on this grid
  int nUnblanksTotal;     //! Total number of cells to unblank across domain
  vector<int> ubCells;    //! Cells from this grid which need to be unblanked
  vector<int> nCellsRecv; //! Number of points incoming from each grid (across
                          //! interComm)
  vector<int>
      nCellsSend; //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>>
      recvCells; //! Cell IDs which will be received from each grid (across
                 //! interComm) (counter to foundPts)

  vector<int> nQptsSend;
  vector<int> nQptsRecv;

  //! Local supermesh of donor elements for each cell needing to be unblanked
  vector<superMesh> donors;

  /* --- Member Functions --- */

  void setup(input *_params, int _nGrids, int _gridID, int _gridRank,
             int _nprocPerGrid, vector<int> &_gridIdList);

  //! Setup the list of overset-boundary flux points to interpolate data to
  void setupOverFacePoints(vector<shared_ptr<overFace>> &overFaces,
                           int nFptsPerFace);

  //! Setup the list of fringe/receptor solution points to interpolate data to
  void setupFringeCellPoints(vector<shared_ptr<ele>> &eles,
                             const unordered_set<int> &fringeCells,
                             const vector<int> &eleMap);

  //! Transfer the exchanged U_in data to the fringe cells
  void transferEleData(vector<shared_ptr<ele>> &eles,
                       const unordered_set<int> &fringeCells,
                       const vector<int> &eleMap);

private:
  /* ---- For Static Cases using Field Interpolation ---- */
  vector<matrix<double>> qpts, qptsD_ref, donorBasis, massMatTDRow, ubLHS;
  vector<vector<int>> targetID, donorID;
  vector<vector<int>> recvInds;

  //! For use with ADT in 2D
  vector<int> eleList;
};
