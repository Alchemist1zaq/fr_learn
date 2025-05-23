
#include "superMesh.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>

#include "funcs.hpp"
#include "points.hpp"

/*
 * Tetrahedron Node-Ordering Conventions:
 *
 *            2                                     2
 *          ,/|`\                                 ,/|`\
 *        ,/  |  `\                             ,/  |  `\
 *      ,/    '.   `\                         ,6    '.   `5
 *    ,/       |     `\                     ,/       8     `\
 *  ,/         |       `\                 ,/         |       `\
 * 0-----------'.--------1               0--------4--'.--------1
 *  `\.         |      ,/                 `\.         |      ,/
 *     `\.      |    ,/                      `\.      |    ,9
 *        `\.   '. ,/                           `7.   '. ,/
 *           `\. |/                                `\. |/
 *              `3                                    `3
 */

superMesh::superMesh() {}

superMesh::superMesh(vector<point> &_target, Array2D<point> &_donors,
                     int _order, int _nDims) {
  setup(_target, _donors, _order, _nDims);
}

superMesh::superMesh(vector<point> &_target, Array2D<point> &_donors,
                     int _order, int _nDims, int _rank, int _ID) {
  rank = _rank;
  ID = _ID;
  setup(_target, _donors, _order, _nDims);
}

superMesh::~superMesh() {}

void superMesh::setup(vector<point> &_target, Array2D<point> &_donors,
                      int _order, int _nDims) {
  target = _target;
  donors = _donors;
  order = _order;
  nDims = _nDims;

  nv_simp = nDims + 1;

  nDonors = donors.getDim0();

  buildSuperMesh();
}

void superMesh::buildSuperMesh(void) {
  if (nDims == 2)
    buildSuperMeshTri();
  else
    buildSuperMeshTet();
}

void superMesh::buildSuperMeshTri(void) {
  // Step 1: Split the donor hexahedrons into triangles to prepare for clipping
  tris.resize(0);
  parents.resize(0);
  for (int i = 0; i < nDonors; i++) {
    auto tmpTris = splitQuadIntoTris(donors.getRow(i));
    tris.insert(tris.end(), tmpTris.begin(), tmpTris.end());
    parents.insert(parents.end(), tmpTris.size(), i);
  }

  // Step 2: Get the clipping planes from the target-cell faces
  normals.resize(4);

  point xc;
  for (auto &pt : target)
    xc += pt;
  xc /= target.size();

  vector<point> facePts = {target[0], target[1]};
  faces.insertRow(facePts);
  normals[0] = getEdgeNormal(facePts, xc);

  facePts = {target[1], target[2]};
  faces.insertRow(facePts);
  normals[1] = getEdgeNormal(facePts, xc);

  facePts = {target[2], target[3]};
  faces.insertRow(facePts);
  normals[2] = getEdgeNormal(facePts, xc);

  facePts = {target[3], target[0]};
  faces.insertRow(facePts);
  normals[3] = getEdgeNormal(facePts, xc);

  // Step 3: Use the faces to clip the tets
  for (uint i = 0; i < faces.getDim0(); i++) {
    vector<triangle> newTris;
    vector<int> newParents;
    for (uint j = 0; j < tris.size(); j++) {
      triangle tri = tris[j];
      auto tmpTris = clipTri(tri, faces.getRow(i), normals[i]);

      newTris.insert(newTris.end(), tmpTris.begin(), tmpTris.end());
      newParents.insert(newParents.end(), tmpTris.size(), parents[j]);
    }
    tris = newTris;
    parents = newParents;
  }
}

void superMesh::buildSuperMeshTet(void) {
  // Step 1: Split the donor hexahedrons into tets to prepare for clipping
  tets.resize(0);
  parents.resize(0);
  for (int i = 0; i < nDonors; i++) {
    auto tmpTets = splitHexIntoTets(donors.getRow(i));
    tets.insert(tets.end(), tmpTets.begin(), tmpTets.end());
    parents.insert(parents.end(), tmpTets.size(), i);
  }

  // Step 2: Get the clipping planes from the target-cell faces
  normals.resize(6); // Consider slightly-more-accurate approach of splitting
                     // faces into tris instead (12 faces)

  point xc;
  for (auto &pt : target)
    xc += pt;
  xc /= target.size();

  vector<point> facePts = {target[0], target[1], target[2], target[3]};
  faces.insertRow(facePts);
  normals[0] = getFaceNormalQuad(facePts, xc);

  facePts = {target[4], target[5], target[6], target[7]};
  faces.insertRow(facePts);
  normals[1] = getFaceNormalQuad(facePts, xc);

  facePts = {target[0], target[3], target[7], target[4]};
  faces.insertRow(facePts);
  normals[2] = getFaceNormalQuad(facePts, xc);

  facePts = {target[2], target[1], target[5], target[6]};
  faces.insertRow(facePts);
  normals[3] = getFaceNormalQuad(facePts, xc);

  facePts = {target[0], target[1], target[5], target[4]};
  faces.insertRow(facePts);
  normals[4] = getFaceNormalQuad(facePts, xc);

  facePts = {target[3], target[2], target[6], target[7]};
  faces.insertRow(facePts);
  normals[5] = getFaceNormalQuad(facePts, xc);

  // Step 3: Use the faces to clip the tets
  for (uint i = 0; i < faces.getDim0(); i++) {
    vector<tetra> newTets;
    vector<int> newParents;
    for (uint j = 0; j < tets.size(); j++) {
      tetra tet = tets[j];
      auto tmpTets = clipTet(tet, faces.getRow(i), normals[i]);
      newTets.insert(newTets.end(), tmpTets.begin(), tmpTets.end());
      newParents.insert(newParents.end(), tmpTets.size(), parents[j]);
    }
    tets = newTets;
    parents = newParents;
  }
}

double superMesh::integrate(const vector<double> &data) {
  if ((int)data.size() != nQpts)
    fatalError(
        "To integrate over supermesh, data must lie at its quadrature nodes.");

  double val = 0;
  for (int i = 0; i < nSimps; i++) {
    for (int j = 0; j < nQpts_simp; j++) {
      val += data[i * nQpts_simp + j] * weights[j] * vol[i];
    }
  }

  return val;
}

vector<double> superMesh::integrate(matrix<double> &data) {
  if ((int)data.getDim0() != nQpts)
    fatalError(
        "To integrate over supermesh, data must lie at its quadrature nodes.");

  int nFields = data.getDim1();

  vector<double> vals(nFields);
  for (int i = 0; i < nSimps; i++)
    for (int j = 0; j < nQpts_simp; j++)
      for (int k = 0; k < nFields; k++)
        vals[k] += data(i * nQpts_simp + j, k) * weights[j] * vol[i];

  return vals;
}

vector<double> superMesh::integrateByDonor(const vector<double> &data) {
  if ((int)data.size() != nQpts)
    fatalError(
        "To integrate over supermesh, data must lie at its quadrature nodes.");

  vector<double> val(nDonors);
  for (int i = 0; i < nSimps; i++) {
    int id = parents[i];
    for (int j = 0; j < nQpts_simp; j++) {
      val[id] += data[i * nQpts_simp + j] * weights[j] * vol[i];
    }
  }

  return val;
}

matrix<double> superMesh::integrateByDonor(matrix<double> &data) {
  if ((int)data.getDim0() != nQpts)
    fatalError(
        "To integrate over supermesh, data must lie at its quadrature nodes.");

  int nFields = data.getDim1();

  matrix<double> val(nDonors, nFields);
  for (int i = 0; i < nSimps; i++) {
    int id = parents[i];
    for (int j = 0; j < nQpts_simp; j++) {
      for (int k = 0; k < nFields; k++) {
        val(id, k) += data(i * nQpts_simp + j, k) * weights[j] * vol[i];
      }
    }
  }

  return val;
}

void superMesh::setupQuadrature(void) {
  if (nDims == 2) {
    getQuadRuleTri(order, qpts, weights);
    /// for (auto &wt:weights) wt /= 2.; // Area of ref. triangle = 1/2
  } else {
    getQuadRuleTet(order, qpts, weights);
    for (auto &wt : weights)
      wt /= 6.; // Volume of ref. tetrahedron = 1/6
  }

  nSimps = (nDims == 2) ? tris.size() : tets.size();
  nQpts_simp = qpts.size();
  nQpts = nQpts_simp * nSimps;

  shapeQpts.setup(nQpts_simp, nv_simp);
  if (nDims == 2) {
    for (int i = 0; i < nQpts_simp; i++)
      shape_tri(qpts[i], shapeQpts[i]);
  } else {
    for (int i = 0; i < nQpts_simp; i++)
      shape_tet(qpts[i], shapeQpts[i]);
  }

  vol.resize(0);

  if (nDims == 2) {
    for (auto &tri : tris) {
      vol.push_back(getAreaTri(tri.nodes));
      tri.qpts.resize(nQpts_simp);
      for (int j = 0; j < nQpts_simp; j++) {
        for (int k = 0; k < nv_simp; k++) {
          tri.qpts[j] += tri.nodes[k] * shapeQpts(j, k);
        }
      }
    }
  } else {
    for (auto &tet : tets) {
      vol.push_back(getVolumeTet(tet.nodes));
      tet.qpts.resize(nQpts_simp);
      for (int j = 0; j < nQpts_simp; j++) {
        for (int k = 0; k < nv_simp; k++) {
          tet.qpts[j] += tet.nodes[k] * shapeQpts(j, k);
        }
      }
    }
  }
}

void superMesh::getQpts(vector<point> &qptPos, vector<int> &qptCell) {
  qptPos.resize(0);
  qptCell.resize(0);
  if (nDims == 2) {
    for (int i = 0; i < nSimps; i++) {
      for (int j = 0; j < nQpts_simp; j++) {
        qptPos.push_back(tris[i].qpts[j]);
        qptCell.push_back(parents[i]);
      }
    }
  } else {
    for (int i = 0; i < nSimps; i++) {
      for (int j = 0; j < nQpts_simp; j++) {
        qptPos.push_back(tets[i].qpts[j]);
        qptCell.push_back(parents[i]);
      }
    }
  }
}

void superMesh::getQpts(matrix<double> &qptPos, vector<int> &qptCell) {
  qptPos.setup(0, 0);
  qptCell.resize(0);
  if (nDims == 2) {
    for (int i = 0; i < nSimps; i++) {
      for (int j = 0; j < nQpts_simp; j++) {
        qptPos.insertRow({tris[i].qpts[j].x, tris[i].qpts[j].y, 0.});
        qptCell.push_back(parents[i]);
      }
    }
  } else {
    for (int i = 0; i < nSimps; i++) {
      for (int j = 0; j < nQpts_simp; j++) {
        qptPos.insertRow(
            {tets[i].qpts[j].x, tets[i].qpts[j].y, tets[i].qpts[j].z});
        qptCell.push_back(parents[i]);
      }
    }
  }
}

vector<double> superMesh::getWeights(void) {
  vector<double> wts = weights;
  return wts;
}

void superMesh::printSuperMesh(int rank, int ID) {
  string filename = "mesh_" + std::to_string((long long)rank) + "_" +
                    std::to_string((long long)ID) + ".csv";
  ofstream mesh(filename.c_str());

  if (nDims == 2) {
    mesh << "Parent,Triangle,x,y" << endl;
    for (int i = 0; i < tris.size(); i++) {
      int id = parents[i];
      point pt = tris[i].nodes[0];
      mesh << id << "," << i << "," << pt.x << "," << pt.y << endl;
      pt = tris[i].nodes[1];
      mesh << id << "," << i << "," << pt.x << "," << pt.y << endl;
      pt = tris[i].nodes[2];
      mesh << id << "," << i << "," << pt.x << "," << pt.y << endl;
    }
  } else {
    mesh << "Parent,Tet,x,y,z" << endl;
    for (int i = 0; i < tets.size(); i++) {
      int id = parents[i];
      point pt = tets[i].nodes[0];
      cout.setf(ios::fixed);
      cout.precision(8);
      mesh << id << ", " << i << ", " << right << setw(14) << pt.x << ","
           << right << setw(14) << pt.y << "," << right << setw(14) << pt.z
           << endl;
      pt = tets[i].nodes[1];
      mesh << id << ", " << i << ", " << right << setw(14) << pt.x << ","
           << right << setw(14) << pt.y << "," << right << setw(14) << pt.z
           << endl;
      pt = tets[i].nodes[2];
      mesh << id << ", " << i << ", " << right << setw(14) << pt.x << ","
           << right << setw(14) << pt.y << "," << right << setw(14) << pt.z
           << endl;
      pt = tets[i].nodes[3];
      mesh << id << ", " << i << ", " << right << setw(14) << pt.x << ","
           << right << setw(14) << pt.y << "," << right << setw(14) << pt.z
           << endl;
    }
  }
  mesh.close();
}

void superMesh::printQuadPoints(vector<double> &vals) {
  string filename = "points_" + std::to_string((long long)rank) + "_" +
                    std::to_string((long long)ID) + ".csv";
  ofstream mesh(filename.c_str());

  if (vals.size() != nQpts) {
    vals.assign(nQpts, 1);
  }

  mesh << "Parent,Triangle,value,x,y" << endl;
  for (int i = 0; i < tris.size(); i++) {
    int id = parents[i];
    for (int j = 0; j < nQpts_simp; j++) {
      point pt = tris[i].qpts[j];
      mesh << id << "," << i << "," << vals[i * nQpts_simp + j] << "," << pt.x
           << "," << pt.y << endl;
    }
  }
  mesh.close();
}

vector<tetra> splitHexIntoTets(const vector<point> &hexNodes) {
  vector<tetra> newTets(5);

  short ind[5][4] = {
      {0, 1, 4, 3}, {2, 1, 6, 3}, {5, 1, 6, 4}, {7, 3, 4, 6}, {1, 3, 6, 4}};

  for (short i = 0; i < 5; i++)
    for (short j = 0; j < 4; j++)
      newTets[i].nodes[j] = hexNodes[ind[i][j]];

  return newTets;
}

vector<tetra> clipTet(tetra &tet, const vector<point> &clipFace, Vec3 &norm) {
  /* --- WARNING: Assuming a linear, planar face --- */

  vector<tetra> outTets;

  // Get face centroid
  point xc;
  for (auto &pt : clipFace)
    xc += pt;
  xc /= clipFace.size();

  set<int> deadPts;
  for (int i = 0; i < 4; i++) {
    // Check each point of tetra to see which must be removed
    Vec3 dx = tet.nodes[i] - xc;
    double dot = dx * norm;
    if (dot > 0) // Point lies on cut-side of clipping plane
      deadPts.insert(i);
  }

  /*
   * Perform the clipping and subdivide the new volume into new tets
   * Only 3 cases in which the clipping can occur
   * New points are created at the intersections of the original tet's edges
   * with the clipping plane: http://geomalgorithms.com/a05-_intersect-1.html
   */
  switch (deadPts.size()) {
  case 0: {
    // No intersection
    outTets.push_back(tet);
    break;
  }

  case 1: {
    // Remove 1 point to get a prism; split prism into 3 new tets
    int kill = *(deadPts.begin()); // The point to remove

    // Get the new points by intersecting the tet's edges with the clipping
    // plane Have to be careful about orientation of final tet
    map<int, array<int, 3>> flipTet;
    flipTet[0] = {{1, 3, 2}};
    flipTet[1] = {{0, 2, 3}};
    flipTet[2] = {{0, 3, 1}};
    flipTet[3] = {{0, 1, 2}};
    array<int, 3> ePts = flipTet[kill];

    // Find the intersection points
    array<point, 3> newPts;
    for (int i = 0; i < 3; i++) {
      Vec3 ab = tet.nodes[ePts[i]] - tet.nodes[kill];
      Vec3 ac = xc - tet.nodes[kill];
      newPts[i] = ab * ((norm * ac) / (norm * ab)) + tet.nodes[kill];
    }

    outTets.resize(3);
    outTets[0].nodes = {{tet.nodes[ePts[0]], tet.nodes[ePts[1]], newPts[0],
                         tet.nodes[ePts[2]]}};
    outTets[1].nodes = {{tet.nodes[ePts[2]], newPts[0], newPts[2], newPts[1]}};
    outTets[2].nodes = {
        {tet.nodes[ePts[1]], tet.nodes[ePts[2]], newPts[1], newPts[0]}};
    break;
  }

  case 2: {
    // Tet cut in half through 4 edges; split into 3 new tets
    // Get the points we're keeping and 'killing'
    array<int, 2> kill, keep;
    int m = 0, n = 0;
    for (int i = 0; i < 4; i++) {
      if (deadPts.count(i)) {
        kill[m] = i;
        m++;
      } else {
        keep[n] = i;
        n++;
      }
    }

    /* Re-orient tet (shuffle nodes) based on kept nodes so that
     * clipping becomes standardized; 'base case' is keeping {0,1}
     * One possible case for each edge being removed */
    map<array<int, 2>, array<int, 4>> flipTet;
    flipTet[{{0, 1}}] = {{0, 1, 2, 3}};
    flipTet[{{0, 2}}] = {{0, 2, 3, 1}};
    flipTet[{{0, 3}}] = {{0, 3, 1, 2}};
    flipTet[{{1, 2}}] = {{1, 2, 0, 3}};
    flipTet[{{1, 3}}] = {{1, 3, 2, 0}};
    flipTet[{{2, 3}}] = {{2, 3, 0, 1}};
    array<int, 4> ind = flipTet[keep];

    // Intersect the plane with the edges to get the new points
    array<point, 4> newPts;
    point a, b;
    Vec3 ab, ac;

    // Edge 0-3
    a = tet.nodes[ind[0]];
    b = tet.nodes[ind[3]];
    ab = b - a;
    ac = xc - a;
    newPts[0] = ab * ((norm * ac) / (norm * ab)) + a;

    // Edge 1-3
    a = tet.nodes[ind[1]];
    b = tet.nodes[ind[3]];
    ab = b - a;
    ac = xc - a;
    newPts[1] = ab * ((norm * ac) / (norm * ab)) + a;

    // Edge 1-2
    a = tet.nodes[ind[1]];
    b = tet.nodes[ind[2]];
    ab = b - a;
    ac = xc - a;
    newPts[2] = ab * ((norm * ac) / (norm * ab)) + a;

    // Edge 0-2
    a = tet.nodes[ind[0]];
    b = tet.nodes[ind[2]];
    ab = b - a;
    ac = xc - a;
    newPts[3] = ab * ((norm * ac) / (norm * ab)) + a;

    // Setup the new tets
    outTets.resize(3);
    outTets[0].nodes = {
        {tet.nodes[ind[1]], newPts[0], newPts[3], tet.nodes[ind[0]]}};
    outTets[1].nodes = {{newPts[0], newPts[3], newPts[1], tet.nodes[ind[1]]}};
    outTets[2].nodes = {{newPts[1], newPts[3], newPts[2], tet.nodes[ind[1]]}};
    break;
  }

  case 3: {
    // The opposite of case 1; new tet is one corner of original tet
    int keep = -1;
    for (int i = 0; i < 4; i++) {
      if (!deadPts.count(i)) {
        keep = i;
        break;
      }
    }

    // Get the new points by intersecting the tet's edges with the clipping
    // plane Have to be careful about orientation of final tet, so map to a
    // 'standard' orientation
    map<int, array<int, 3>> flipTet;
    flipTet[0] = {{1, 3, 2}};
    flipTet[1] = {{0, 2, 3}};
    flipTet[2] = {{0, 3, 1}};
    flipTet[3] = {{0, 1, 2}};
    array<int, 3> ePts = flipTet[keep];

    // Setup outgoing tet; node 3 is the 'kept' node
    // Find the intersection points
    outTets.resize(1);
    outTets[0].nodes[3] = tet.nodes[keep];
    for (int i = 0; i < 3; i++) {
      Vec3 ab = tet.nodes[ePts[i]] - tet.nodes[keep];
      Vec3 ac = xc - tet.nodes[keep];
      outTets[0].nodes[i] = ab * ((norm * ac) / (norm * ab)) + tet.nodes[keep];
    }
    break;
  }

  case 4: {
    // Entire tet is beyond clipping face
    break;
  }
  }

  return outTets;
}

vector<triangle> clipTri(triangle &tri, const vector<point> &clipEdge,
                         Vec3 &norm) {
  /* --- WARNING: Assuming a linear edge --- */

  vector<triangle> outTris;

  // Get face centroid
  point xc = clipEdge[0];
  xc += clipEdge[1];
  xc /= 2.;

  set<int> deadPts;
  for (int i = 0; i < 3; i++) {
    // Check each point of triangle to see which must be removed
    Vec3 dx = tri.nodes[i] - xc;
    double dot = dx * norm;
    if (dot > 0) // Point lies on cut-side of clipping plane
      deadPts.insert(i);
  }

  /*
   * Perform the clipping and subdivide the new volume into new tris
   */
  switch (deadPts.size()) {
  case 0: {
    // No intersection.
    outTris.push_back(tri);
    break;
  }
  case 1: {
    // Removing one corner of tri
    int kill = *(deadPts.begin()); // The point to remove

    // Get the points being kept; map to a 'standard' triangle
    map<int, array<int, 2>> flipTri;
    flipTri[0] = {{1, 2}};
    flipTri[1] = {{2, 0}};
    flipTri[2] = {{0, 1}};
    array<int, 2> ePts = flipTri[kill];

    // Find the cutting-plane intersection points
    Vec3 ab = tri.nodes[ePts[0]] - tri.nodes[kill];
    Vec3 ac = xc - tri.nodes[kill];
    point newPt1 = ab * (norm * ac) / (norm * ab) + tri.nodes[kill];

    ab = tri.nodes[ePts[1]] - tri.nodes[kill];
    point newPt2 = ab * (norm * ac) / (norm * ab) + tri.nodes[kill];

    outTris.resize(2);
    outTris[0].nodes = {{tri.nodes[ePts[0]], tri.nodes[ePts[1]], newPt1}};
    outTris[1].nodes = {{tri.nodes[ePts[1]], newPt2, newPt1}};
    break;
  }
  case 2: {
    // Keeping one corner of tri
    int keep = -1;
    for (int i = 0; i < 3; i++) {
      if (!deadPts.count(i)) {
        keep = i;
        break;
      }
    }

    map<int, array<int, 2>> flipTri;
    flipTri[0] = {{1, 2}};
    flipTri[1] = {{2, 0}};
    flipTri[2] = {{0, 1}};
    array<int, 2> ePts = flipTri[keep];

    // Setup outgoing tri; node 2 is the 'kept' node
    // Find the intersection points
    outTris.resize(1);
    outTris[0].nodes[2] = tri.nodes[keep];
    for (int i = 0; i < 2; i++) {
      Vec3 ab = tri.nodes[ePts[i]] - tri.nodes[keep];
      Vec3 ac = xc - tri.nodes[keep];
      outTris[0].nodes[i] = ab * ((norm * ac) / (norm * ab)) + tri.nodes[keep];
    }
    break;
  }
  case 3: {
    // Entire tri is beyond clipping face
    break;
  }
  }

  return outTris;
}

vector<triangle> splitQuadIntoTris(const vector<point> &quadNodes) {
  vector<triangle> newTris(2);

  newTris[0].nodes = {{quadNodes[0], quadNodes[1], quadNodes[3]}};
  newTris[1].nodes = {{quadNodes[1], quadNodes[2], quadNodes[3]}};

  return newTris;
}

double getAreaTri(std::array<point, 3> &nodes) {
  Vec3 dx1 = nodes[1] - nodes[0];
  Vec3 dx2 = nodes[2] - nodes[0];
  Vec3 cross = dx1.cross(dx2);

  return 1. / 2. * std::abs(cross.z);
}

double getVolumeTet(std::array<point, 4> &nodes) {
  Vec3 dx1 = nodes[1] - nodes[0];
  Vec3 dx2 = nodes[2] - nodes[0];
  Vec3 dx3 = nodes[3] - nodes[0];
  Vec3 cross = dx1.cross(dx2);

  return 1. / 6. * std::abs(dx3 * cross);
}
