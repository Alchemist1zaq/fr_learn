
#pragma once

#include <array>
#include <vector>

#include "global.hpp"

#include "matrix.hpp"

// Possibility: use 'class Simplex'; 'class tetra: public Simplex';
// 'vector<shared_ptr<Simplex>>' ...

struct tetra {
  array<point, 4> nodes; //! Positions of nodes in tet
  vector<point> qpts;    //! Physical positions of quadrature points in tet
  // vector<double> weights;  //! Weights associated with each quadrature point
  int donorID; //! Donor-grid cell ID to which tetra belongs
};

struct triangle {
  array<point, 3> nodes; //! Positions of nodes in tri
  vector<point> qpts;    //! Physical positions of quadrature points in tri
  int donorID;           //! Donor-grid cell ID to which tetra belongs
};

class superMesh {
public:
  /* ---- Default Functions ---- */

  superMesh();
  ~superMesh();

  superMesh(vector<point> &_target, Array2D<point> &_donors, int _order,
            int _nDims);
  superMesh(vector<point> &_target, Array2D<point> &_donors, int _order,
            int _nDims, int _rank, int _ID);
  /* ---- Member Variables ---- */

  vector<point> target;  //! Target cell's node positions for which to create
                         //! local supermesh
  Array2D<point> donors; //! Node positions of cells from donor grid which
                         //! overlap target cell

  int nSimps;     //! Total number of simplices comprising the supermesh
  int nv_simp;    //! Number of vertices in simplex (3 or 4)
  int nQpts;      //! Total number of quadrature points in the whole supermesh
  int order;      //! Order of quadrature rule to use
  int nQpts_simp; //! Number of quadrature points per tet (based on order)
  int nDonors;    //! Number of donor cells
  int nDims;      //! Dimension of the problem: 2D or 3D

  Array2D<point>
      faces; //! Face points of target cell for use as clipping planes
  vector<Vec3> normals; //! Outward face normals for target cell (for clipping)

  vector<point>
      qpts; //! Locations of quadrature points in reference tetrahedron
  vector<double> weights;   //! Quadrature weights
  matrix<double> shapeQpts; //! Values of tetrahedron nodal shape basis at
                            //! quadrature points [nQpts x 4]

  vector<tetra> tets;    //! Tetrahedrons comprising the supermesh [for 3D]
  vector<triangle> tris; //! triangles comprising the supermesh [for 2D]
  vector<int>
      parents; //! Parent donor-cell ID for each tet [range 0:(nDonors-1)]
  vector<double> vol; //! Volume of each simplex

  /* ---- Member Functions ---- */

  void setup(vector<point> &_target, Array2D<point> &_donors, int _order,
             int nDims);

  //! Using given grids and target cell, build the local supermesh
  void buildSuperMesh(void);

  //! Integrate a quantity over the superMesh (given at quadrature points of
  //! supermesh tets)
  double integrate(const vector<double> &data);

  //! Integrate a series of quantities over the superMesh (given at quadrature
  //! points of supermesh tets)
  vector<double> integrate(matrix<double> &data);

  //! Integrate the given data over each donor individually
  vector<double> integrateByDonor(const vector<double> &data);

  //! Integrate each field of the given data over each donor individually
  matrix<double> integrateByDonor(matrix<double> &data);

  /*!
   * Returns the physical(?) positions of the quadrature points, along
   * with which donor-grid cell they lie within
   */
  void getQpts(vector<point> &qptPos, vector<int> &qptCell);
  void getQpts(matrix<double> &qptPos, vector<int> &qptCell);

  vector<double> getWeights(void);

  int getNQpts(void) { return nQpts; }

  //! Get the quadrature point locations and weights, and find physical
  //! positions of all qpts
  void setupQuadrature(void);

  //! Print the simplices of the superMesh to a CSV file
  void printSuperMesh(int rank, int ID);

  //! Print out the x,y positions on all quadrature points, along with values at
  //! each point
  void printQuadPoints(vector<double> &vals);

  int ID, rank;

private:
  void buildSuperMeshTri(void);
  void buildSuperMeshTet(void);
};

/* --- Extra Helper Functions --- */

//! Subdivide the given hexahedron into 5 tetrahedrons
vector<tetra> splitHexIntoTets(const vector<point> &hexNodes);

//! Subdivide the given quadrilateral into 2 triangles
vector<triangle> splitQuadIntoTris(const vector<point> &quadNodes);

//! Use the given face and outward normal to clip the given tet and return the
//! new set of tets
vector<tetra> clipTet(tetra &tet, const vector<point> &clipFace, Vec3 &norm);

//! Use the given face and outward normal to clip the given triangle and return
//! the new set of tris
vector<triangle> clipTri(triangle &tri, const vector<point> &clipFace,
                         Vec3 &norm);

double getAreaTri(std::array<point, 3> &nodes);

double getVolumeTet(std::array<point, 4> &nodes);
