
#pragma once

#include <array>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "global.hpp"

#include "ele.hpp"
#include "face.hpp"
#include "input.hpp"
#include "solver.hpp"

#define NORMAL 1
#define HOLE 0
#define FRINGE -1
#define FIELD_HOLE -2

class geo {
public:
  geo();

  ~geo();

  /* === Primay setup routines === */

  //! Setup the geomery using input parameters
  void setup(input *params, bool HMG = false);

  //! Multigrid-specific setup function [from mesh-refinement method]
  void setup_hmg(input *params, int _gridID, int _gridRank, int _nProcGrid,
                 const vector<int> &_gridIdList = {0});

  //! Take the basic connectivity data and generate the rest
  void processConnectivity();

  //! Create the elements and faces needed for the simulation
  void setupElesFaces(input *params, vector<shared_ptr<ele>> &eles,
                      vector<shared_ptr<face>> &faces);

  //! Read essential connectivity from a Gmsh mesh file
  void readGmsh(string fileName);

  int nDims, nFields;
  int nEles, number_of_vertexs, nEdges, nFaces, nIntFaces, nBndFaces, nMpiFaces;
  int nBounds; //! Number of boundaries
  int meshType;
  int nNodesPerCell;

  // Basic [essential] Connectivity Data
  matrix<int> c2v;
  matrix<double>
      xv; //! Current physical position of vertices [static or moving grids]

  // Basic Moving-Grid Variables
  vector<point> xv_new; //! Physical position of vertices for next time step
                        //! [moving grids]
  vector<point> xv0;    //! Initial position of vertices [moving grids]
  matrix<double> rv0;   //! Initial posiiton of vertices in polar coords [moving
                        //! grids, vibrating sphere test case]
  matrix<double> gridVel; //! Grid velocity of vertices

  point minPt; //! Centroid of all vertices on grid partition
  point maxPt; //! Overall x,y,z extents (max-min) of grid partition

  // Additional Connectivity Data
  matrix<int> c2e, c2b, e2c, e2v, v2e, v2v, v2c;
  matrix<int> c2f, f2v, f2c, c2c, c2ac;
  vector<int> v2nv, v2nc, c2nv, c2nf, f2nv, ctype;
  vector<int> intFaces, bndFaces, mpiCells;
  unordered_set<int> overFaces,
      overCells;          //! List of all faces / cells which have an
                          //! overset-boundary-condition face
  vector<int> bcList;     //! List of boundary conditions for each boundary
  vector<string> bcNames; //! List of boundaries given in mesh file
  vector<int> bcType;     //! Boundary condition for each boundary face
  matrix<int> bndPts;     //! List of node IDs on each boundary

  vector<int> bcID;
  int number_of_boundaries;
  vector<int> nBndPts; //! Number of points on each boudary
  vector<matrix<int>>
      bcFaces; //! List of nodes on each face (edge) for each boundary condition
  vector<int> nFacesPerBnd; //! List of # of faces on each boundary
  vector<int> procR;        //! What processor lies to the 'right' of this face

  vector<int>
      gIC_R; //! The global cell ID of the right cell on the opposite processor
  vector<int> mpiLocF;   //! Element-local face ID of MPI Face in left cell
  vector<int> mpiLocF_R; //! Element-local face ID of MPI Face in right cell
  vector<int>
      mpiPeriodic;      //! Flag for whether an MPI face is also a periodic face
  vector<int> faceType; //! Type for each face: hole, internal, boundary, MPI,
                        //! overset [-1,0,1,2,3]

  /* --- Overset-Related Variables --- */
  int nGrids;    //! Number of distinct overset grids
  int nProcGrid; //! Number of MPI processes assigned to current (overset) grid
                 //! block
  int gridID;    //! Which (overset) grid block is this process handling
  int gridRank;  //! MPI rank of process *within* the grid block [0 to
                 //! nprocPerGrid-1]
  int rank;
  int nproc;
  // vector<int> nProcsGrid; //! Number of processes for each (overset) grid
  // block
  vector<int> gridIdList; //! gridID for each MPI rank
  vector<int> iblank; //! Output of TIOGA: flag for whether vertex is normal,
                      //! blanked, or receptor
  vector<int> iblankCell; //! Output? of TIOGA: flag for whether cell is normal,
                          //! blanked, or receptor
  vector<int>
      iblankFace;    //! Flag for whether a face is normal, blanked, or receptor
  vector<int> iwall; //! List of nodes on wall boundaries
  vector<int> iover; //! List of nodes on overset boundaries
  vector<int>
      nodeType; //! For each node: normal interior, normal boundary, or overset

  matrix<int>
      wallFaceNodes; //! For 2D: All the wall-boundary faces for hole cutting

  set<int> mpiNodes; //! Set of all nodes which lie on an MPI boundary

  vector<int> eleMap;  //! For overset meshes where some cells are blanked, map
                       //! from 'ic' to 'eles' index
  vector<int> faceMap; //! For overset meshes where some faces are blanked, map
                       //! from 'ff' to faceType-vector index
  vector<int> currFaceType; //! Current face class type for each face in mesh
                            //! [internal, boundary, mpi, overset, hole]

  vector<int> epart; //! Cell partition data

  /* --- Moving-Overset-Grid-Related Variables --- */
  unordered_set<int>
      holeCells; //! List of cells in mesh which are currently blanked
  unordered_set<int>
      holeFaces; //! List of faces in mesh which are currently blanked
  unordered_set<int> unblankCells; //! List of non-existing cells which, due to
                                   //! motion, must be un-blanked
  unordered_set<int> blankCells;   //! List of existing cells which, due to
                                   //! motion, must be blanked
  unordered_set<int> fringeCells;  //! For field-fill (non-boundary) overset
                                   //! method, fringe/receptor cell list

  int *nodesPerCell; //! Pointer for Tioga to know # of nodes for each element
                     //! type
  array<int *, 1> conn; //! Pointer to c2v for each element type [but only 1, so
                        //! will be size(1)]
  matrix<double> eleBBox;

  void refineGrid2D(geo &outGrid, int nLevels, int shapeOrder);

private:
  input *params;

  /* --- MPI-Related Varialbes (global vs. local data) --- */
  matrix<int> c2v_g;                   //! Global element connectivity
  matrix<double> xv_g;                 //! Global mesh node locations
  vector<int> ic2icg;                  //! Local cell to global cell index
  vector<int> iv2ivg;                  //! Local vertex to global vertex index
  vector<int> ctype_g, c2ne_g, c2nv_g; //! Global element info
  matrix<int> bndPts_g;                //! Global lists of points on boundaries
  vector<int> nBndPts_g; //! Global number of points on each boundary
  map<int, int> bcIdMap; //! Map from Gmsh boundary ID to Flurry BC ID
  int nEles_g, number_of_vertexs_g;

  void processConn3D(void);
  void processConnExtra(void);

  //! Compare the orientation (rotation in ref. space) betwen the local faces of
  //! 2 elements
  int compareOrientation(int ic1, int ic2, int f1, int f2);

  //! Compare two faces [lists of nodes] to see if they match [used for MPI]
  bool compareFaces(vector<int> &face1, vector<int> &face2);
};
