#include "geo.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_set>

#include "boundFace.hpp"
#include "face.hpp"
#include "intFace.hpp"

geo::geo() { nodesPerCell = NULL; }

geo::~geo() {
  if (nodesPerCell != NULL) {
    delete[] nodesPerCell;
  }
}

void geo::setup(input *params, bool HMG) {
  this->params = params;

  nDims = params->nDims;
  nFields = params->nFields;
  meshType = params->meshType;
  gridID = 0;
  gridRank = params->rank;
  nProcGrid = params->nproc;
  rank = params->rank;
  nproc = params->nproc;
  if (nDims == 2) {
    fatalError("no 2 dim");
  }
  readGmsh(params->meshFileName);

  nNodesPerCell = getMax(c2nv);

  if (HMG) {
    processConn3D();
  } else {
    processConnectivity();
  }
}

void geo::setup_hmg(input *params, int _gridID, int _gridRank, int _nProcGrid,
                    const vector<int> &_gridIdList) {
  this->params = params;

  nDims = params->nDims;
  nFields = params->nFields;
  meshType = params->meshType;
  rank = params->rank;
  nproc = params->nproc;
  nGrids = params->nGrids;

  gridID = _gridID;
  gridRank = _gridRank;
  nProcGrid = _nProcGrid;
  gridIdList = _gridIdList;

  processConnectivity();
}

void geo::processConnectivity() {
  cout << "Geo: Processing element connectivity" << endl;

  processConn3D();
  processConnExtra();
}

void geo::processConn3D(void) {
  /* --- Setup Single List of All Faces (sorted vertex lists) --- */

  matrix<int> f2v1, e2v1;
  vector<int> f2nv1;

  // Handy map to store local face-vertex lists for each ele type
  map<int, matrix<int>> ct2fv;
  map<int, vector<int>> ct2fnv;
  // --- FIX ORDERING FOR FUTURE USE ---
  ct2fv[HEX].insertRow(vector<int>{0, 1, 2, 3}); // Bottom (zmin)
  ct2fv[HEX].insertRow(vector<int>{4, 5, 6, 7}); // Top    (zmax)
  ct2fv[HEX].insertRow(vector<int>{3, 0, 4, 7}); // Left   (xmin)
  ct2fv[HEX].insertRow(vector<int>{2, 1, 5, 6}); // Right  (xmax)
  ct2fv[HEX].insertRow(vector<int>{1, 0, 4, 5}); // Front  (ymin)
  ct2fv[HEX].insertRow(vector<int>{3, 2, 6, 7}); // Back   (ymax)
  ct2fnv[HEX] = {4, 4, 4, 4, 4, 4};

  for (int e = 0; e < nEles; e++) {
    for (int f = 0; f < c2nf[e]; f++) {
      // Get local vertex list for face
      auto iface = ct2fv[ctype[e]].getRow(f);

      // Get global vertex list for face
      vector<int> facev(ct2fnv[ctype[e]][f]);
      for (int i = 0; i < ct2fnv[ctype[e]][f]; i++) {
        facev[i] = c2v(e, iface[i]);
        if (i > 0 && facev[i] == facev[i - 1])
          facev[i] = -1;
      }

      // Sort the vertices for easier comparison later
      std::sort(facev.begin(), facev.end());
      f2v1.insertRowUnsized(facev);
      f2nv1.push_back(ct2fnv[ctype[e]][f]);

      e2v1.insertRow({facev[0], facev[1]});
      e2v1.insertRow({facev[1], facev[2]});
      e2v1.insertRow({facev[2], facev[3]});
      e2v1.insertRow({facev[0], facev[3]});
    }
  }

  /* --- Get a unique list of faces --- */

  // NOTE: Could setup f2c here, but I already have another algorithm
  // implemented later

  // iE is of length [original f2v1] with range [final f2v]
  // The number of times a face appears in iF is equal to
  // the number of cells that face touches
  vector<int> iF, iE;
  f2v1.unique(f2v, iF);
  e2v1.unique(e2v, iE);
  nFaces = f2v.getDim0();
  nEdges = e2v.getDim0();

  f2nv.resize(nFaces);
  for (uint i = 0; i < f2nv1.size(); i++)
    f2nv[iF[i]] = f2nv1[i];

  /* --- Generate Internal and Boundary Face Lists --- */

  // Flag for whether global face ID corresponds to interior or boundary face
  // (note that, at this stage, MPI faces will be considered boundary faces)
  faceType.assign(nFaces, INTERNAL);

  nIntFaces = 0;
  nBndFaces = 0;
  nMpiFaces = 0;

  for (uint i = 0; i < iF.size(); i++) {
    if (iF[i] != -1) {
      auto ff = findEq(iF, iF[i]);
      if (ff.size() > 2) {
        stringstream ss;
        ss << i;
        string errMsg = "More than 2 cells for face " + ss.str();
        fatalError(errMsg.c_str());
      } else if (ff.size() == 2) {
        // Internal Edge which has not yet been added
        intFaces.push_back(iF[i]);
        nIntFaces++;
      } else if (ff.size() == 1) {
        // Boundary or MPI Edge
        bndFaces.push_back(iF[i]);
        faceType[iF[i]] = BOUNDARY;
        nBndFaces++;
      }

      // Mark edges as completed
      vecAssign(iF, ff, -1);
    }
  }

  /* --- Match Boundary Faces to Boundary Conditions --- */

  bcFaces.resize(nBounds);
  bcType.assign(nBndFaces, NONE);
  bcID.resize(nBndFaces);
  for (int i = 0; i < nBndFaces; i++) {
    for (int bnd = 0; bnd < nBounds; bnd++) {
      bool isOnBound = true;
      for (int j = 0; j < f2nv[bndFaces[i]]; j++) {
        if (findFirst(bndPts[bnd], f2v(bndFaces[i], j), bndPts.dims[1]) == -1) {
          isOnBound = false;
          break;
        }
      }

      if (isOnBound) {
        // cout << "bndFace matched to bc " << bcList[bnd] << endl;
        //  The edge lies on this boundary
        bcType[i] = bcList[bnd];
        bcID[i] = bnd;
        bcFaces[bnd].insertRow(f2v[bndFaces[i]], INSERT_AT_END, f2v.dims[1]);
        break;
      }
    }
  }

  /* --- Setup Cell-To-Face, Face-To-Cell --- */

  c2f.setup(nEles, getMax(c2nf));
  c2b.setup(nEles, getMax(c2nf));
  c2c.setup(nEles, getMax(c2nf));
  c2f.initializeToValue(-1);
  c2b.initializeToZero();
  c2c.initializeToValue(-1);
  f2c.setup(nFaces, 2);
  f2c.initializeToValue(-1);

  for (int ic = 0; ic < nEles; ic++) {
    for (int j = 0; j < c2nf[ic]; j++) {
      // Get local vertex list for face
      auto iface = ct2fv[ctype[ic]].getRow(j);

      // Get global vertex list for face
      int fnv = ct2fnv[ctype[ic]][j];
      vector<int> facev(fnv);
      for (int i = 0; i < fnv; i++)
        facev[i] = c2v(ic, iface[i]);

      // Sort the vertices for easier comparison
      std::sort(facev.begin(), facev.end());

      bool found = false;

      // Check if face is actually collapsed (nonexistant) (all nodes identical)
      bool collapsed = true;
      for (int i = 1; i < fnv; i++)
        collapsed = (collapsed && (facev[0] == facev[i]));

      if (collapsed)
        continue;

      for (int f = 0; f < nFaces; f++) {
        if (std::equal(f2v[f], f2v[f] + fnv, facev.begin())) {
          found = true;
          c2f(ic, j) = f;
          break;
        }
      }

      if (!found)
        fatalError("Unable to match cell face to global face list!");

      int ff = c2f(ic, j);

      // Find ID of face within type-specific array
      if (faceType[ff] > 0)
        c2b(ic, j) = 1;
      else
        c2b(ic, j) = 0;

      if (f2c(ff, 0) == -1) {
        // No cell yet assigned to edge; put on left
        f2c(ff, 0) = ic;
      } else {
        // Put cell on right
        f2c(ff, 1) = ic;
        // Update c2c for both cells
        int ic2 = f2c(ff, 0);
        if (ic2 != ic) {
          vector<int> cellFaces(c2f[ic2], c2f[ic2] + c2nf[ic2]);
          int fid2 = findFirst(cellFaces, ff);
          c2c(ic, j) = ic2;
          c2c(ic2, fid2) = ic;
        }
      }
    }
  }
}

void geo::processConnExtra(void) {
  int maxNC;
  maxNC = 26;

  getBoundingBox(xv, minPt, maxPt);

  // Get vertex to vertex/edge connectivity
  vector<set<int>> v2v_tmp(number_of_vertexs);
  vector<set<int>> v2e_tmp(number_of_vertexs);
  if (nDims == 2) {
    // Faces are the edges
    for (int ie = 0; ie < nFaces; ie++) {
      int iv1 = f2v(ie, 0);
      int iv2 = f2v(ie, 1);
      v2v_tmp[iv1].insert(iv2);
      v2v_tmp[iv2].insert(iv1);
      v2e_tmp[iv1].insert(ie);
      v2e_tmp[iv2].insert(ie);
    }
  } else {
    // Edges were processed separately for 3D
    for (int ie = 0; ie < nEdges; ie++) {
      int iv1 = e2v(ie, 0);
      int iv2 = e2v(ie, 1);
      v2v_tmp[iv1].insert(iv2);
      v2v_tmp[iv2].insert(iv1);
      v2e_tmp[iv1].insert(ie);
      v2e_tmp[iv2].insert(ie);
    }
  }

  v2nv.resize(number_of_vertexs);
  for (int iv = 0; iv < number_of_vertexs; iv++) {
    v2nv[iv] = v2v_tmp[iv].size();
  }

  v2v.setup(number_of_vertexs, getMax(v2nv));
  v2e.setup(number_of_vertexs, getMax(v2nv));
  for (int iv = 0; iv < number_of_vertexs; iv++) {
    int j = 0;
    for (auto &iv2 : v2v_tmp[iv]) {
      v2v(iv, j) = iv2;
      j++;
    }

    j = 0;
    for (auto &ie : v2e_tmp[iv]) {
      v2e(iv, j) = ie;
      j++;
    }
  }
}

void geo::setupElesFaces(input *params, vector<shared_ptr<ele>> &eles,
                         vector<shared_ptr<face>> &faces) {
  if (nEles <= 0)
    fatalError("Cannot setup elements array - nEles = 0");

  if (params->rank == 0) {
    cout << "Geo: Setting up elements" << endl;
  }

  eles.resize(0);
  faces.resize(0);

  // Setup the elements

  eleMap.assign(nEles, -1);

  int nc = 0;
  for (int ic = 0; ic < nEles; ic++) {
    // Skip any hole cells
    if (meshType == OVERSET_MESH && iblankCell[ic] == HOLE)
      continue;

    shared_ptr<ele> e = make_shared<ele>();
    e->ID = ic;
    e->sID = nc;
    if (nProcGrid > 1)
      e->IDg = ic2icg[ic];
    else
      e->IDg = ic;
    e->eType = ctype[ic];
    e->nNodes = c2nv[ic];
    if (nDims == 2)
      e->nMpts = 4;
    else
      e->nMpts = 8;

    eles.push_back(e);
    eleMap[ic] = nc;
    nc++;
  }

  /* --- Setup the faces --- */
  vector<int> cellFaces;

  if (params->rank == 0) {
    cout << "Geo: Setting up internal faces" << endl;
  }

  faceMap.assign(nFaces, -1);
  currFaceType.assign(nFaces, HOLE_FACE);

  // Internal Faces
  for (auto &ff : intFaces) {
    // Skip any hole faces
    if (meshType == OVERSET_MESH && iblankFace[ff] != NORMAL)
      continue;

    shared_ptr<face> iface = make_shared<intFace>();

    int ic1 = f2c(ff, 0);
    // Find local face ID of global face within first element [on left]
    cellFaces.assign(c2f[ic1], c2f[ic1] + c2nf[ic1]);
    int fid1 = findFirst(cellFaces, ff);
    if (f2c(ff, 1) == -1) {
      fatalError("Interior face does not have a right element assigned.");
    } else {
      int ic2 = f2c(ff, 1);
      cellFaces.assign(c2f[ic2], c2f[ic2] + c2nf[ic2]); // List of cell's faces
      int fid2 = findFirst(cellFaces, ff); // Which one is this face
      int relRot = compareOrientation(ic1, fid1, ic2, fid2);
      struct faceInfo info;
      info.IDR = fid2;
      info.relRot = relRot;
      ic1 = eleMap[ic1];
      ic2 = eleMap[ic2];
      if (ic2 == -1)
        fatalError("Internal face has right cell blanked.");
      iface->initialize(eles[ic1], eles[ic2], ff, fid1, info, params);
    }

    faces.push_back(iface);
    faceMap[ff] = faces.size() - 1;
    currFaceType[ff] = INTERNAL;
  }
  // New # of internal faces (after excluding blanked faces)
  nIntFaces = faces.size();

  if (params->rank == 0) {
    cout << "Geo: Setting up boundary faces" << endl;
  }

  // Boundary Faces
  for (int i = 0; i < nBndFaces; i++) {
    // Find global face ID of current boundary face
    int ff = bndFaces[i];

    if ((meshType == OVERSET_MESH && iblankFace[ff] != NORMAL) ||
        bcType[i] == OVERSET)
      continue;

    shared_ptr<face> bface = make_shared<boundFace>();

    int ic = f2c(ff, 0);
    // Find local face ID of global face within element
    cellFaces.assign(c2f[ic], c2f[ic] + c2nf[ic]);
    int fid1 = findFirst(cellFaces, ff);
    if (f2c(ff, 1) != -1) {
      fatalError("Boundary face has a right element assigned.");
    } else {
      struct faceInfo info;
      info.bcType = bcType[i];
      info.isBnd = 1;
      ic = eleMap[ic];
      shared_ptr<ele>
          nullEle; // Since just giving the funciton 'NULL' isn't possible
      bface->initialize(eles[ic], nullEle, ff, fid1, info, params);
    }

    faces.push_back(bface);
    faceMap[ff] = faces.size() - 1;
    currFaceType[ff] = BOUNDARY;
  }
  // New # of boundary faces (after excluding blanked faces)
  nBndFaces = faces.size() - nIntFaces;
}

void geo::readGmsh(string fileName) {
  ifstream meshFile;
  string str;

  cout << "Geo: Reading mesh file " << fileName << endl;

  meshFile.open(fileName.c_str());
  if (!meshFile.is_open())
    fatalError("Unable to open mesh file.");

  // Move cursor to $PhysicalNames
  while (1) {
    getline(meshFile, str);
    if (str.find("$PhysicalNames") != string::npos) {
      break;
    }
    if (meshFile.eof()) {
      fatalError("$PhysicalNames tag not found in Gmsh file!");
    }
  }

  // Read number of boundaries and fields defined
  meshFile >> number_of_boundaries;
  getline(meshFile, str); // clear rest of line

  nBounds = 0;
  for (int i = 0; i < number_of_boundaries; i++) {
    string bcStr, bcName;
    stringstream ss;
    int bcdim, bcid;

    getline(meshFile, str);
    ss << str;
    ss >> bcdim >> bcid >> bcStr;

    // Remove quotation marks ("") from around boundary condition
    size_t ind = bcStr.find("\"");
    while (ind != string::npos) {
      bcStr.erase(ind, 1);
      ind = bcStr.find("\"");
    }
    bcName = bcStr;

    // Convert to lowercase to match boundary condition strings
    std::transform(bcStr.begin(), bcStr.end(), bcStr.begin(), ::tolower);

    // First, map mesh boundary to boundary condition in input file
    if (!params->meshBounds.count(bcStr)) {
      string errS = "Unrecognized mesh boundary: \"" + bcStr + "\"\n";
      errS += "Boundary names in input file must match those in mesh file.";
      fatalError(errS.c_str());
    }

    // Map the Gmsh PhysicalName to the input-file-specified Flurry boundary
    // condition
    bcStr = params->meshBounds[bcStr];

    // Next, check that the requested boundary condition exists
    if (!bcStr2Num.count(bcStr)) {
      string errS = "Unrecognized boundary condition: \"" + bcStr + "\"";
      fatalError(errS.c_str());
    }

    if (bcStr.compare("fluid") == 0) {
      nDims = bcdim;
      params->nDims = bcdim;
      bcIdMap[bcid] = -1;
    } else {
      bcList.push_back(bcStr2Num[bcStr]);
      bcNames.push_back(bcName);
      bcIdMap[bcid] = nBounds; // Map Gmsh bcid to Flurry bound index
      nBounds++;
    }
  }

  // --- Read Mesh Vertex Locations ---
  // Move cursor to $Nodes
  meshFile.clear();
  meshFile.seekg(0, ios::beg);
  while (1) {
    getline(meshFile, str);
    if (str.find("$Nodes") != string::npos) {
      {
        break;
      }
    }
    if (meshFile.eof()) {
      fatalError("$Nodes tag not found in Gmsh file!");
    }
  }

  meshFile >> number_of_vertexs;
  xv.setup(number_of_vertexs, nDims);
  getline(meshFile, str); // Clear end of line, just in case

  for (int i = 0; i < number_of_vertexs; i++) {
    uint32_t iv;
    meshFile >> iv >> xv(i, 0) >> xv(i, 1) >> xv(i, 2);
    getline(meshFile, str);
  }

  // --- Read Element Connectivity ---
  // Move cursor to $Elements
  meshFile.clear();
  meshFile.seekg(0, ios::beg);
  while (1) {
    getline(meshFile, str);
    if (str.find("$Elements") != string::npos) {
      break;
    }
    if (meshFile.eof()) {
      fatalError("$Elements tag not found in Gmsh file!");
    }
  }

  int number_of_elements;
  vector<int> c2v_tmp(27, 0); // Maximum number of nodes/element possible
  vector<set<int>> boundPoints(nBounds);

  nBndPts.resize(nBounds);

  // Read total number of interior + boundary elements
  meshFile >> number_of_elements;
  getline(meshFile, str); // Clear end of line, just in case

  // For Gmsh node ordering, see:
  // http://geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
  int ic = 0;
  for (int k = 0; k < number_of_elements; k++) {
    int id, eType, nTags, bcid, tmp;
    meshFile >> id >> eType >> nTags;
    meshFile >> bcid;
    int gmshID = bcid;
    bcid = bcIdMap[bcid];
    for (int tag = 0; tag < nTags - 1; tag++)
      meshFile >> tmp;

    if (bcid == -1) {
      // NOTE: Currently, only Hexahedron are supported
      switch (eType) {
      case 5:
        // Linear Hexahedron
        c2nv.push_back(8);
        c2nf.push_back(6);
        ctype.push_back(HEX);
        for (int i = 0; i < 8; i++)
          meshFile >> c2v_tmp[i];
        break;

      case 17:
        // Quadratic (20-Node Serendipity) Hexahedron
        c2nv.push_back(20);
        c2nf.push_back(6);
        ctype.push_back(HEX);
        // Corner Nodes
        meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >>
            c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
        // Edge Nodes
        meshFile >> c2v_tmp[8] >> c2v_tmp[11] >> c2v_tmp[12] >> c2v_tmp[9] >>
            c2v_tmp[13] >> c2v_tmp[10];
        meshFile >> c2v_tmp[14] >> c2v_tmp[15] >> c2v_tmp[16] >> c2v_tmp[19] >>
            c2v_tmp[17] >> c2v_tmp[18];
        break;

      case 12:
        // Quadratic (27-Node Lagrange) Hexahedron
        c2nv.push_back(27);
        c2nf.push_back(6);
        ctype.push_back(HEX);
        c2v_tmp.resize(27);
        for (int i = 0; i < c2nv.back(); i++)
          meshFile >> c2v_tmp[i];
        break;

      case 92:
        // Cubic Hexahedron
        c2nv.push_back(64);
        c2nf.push_back(6);
        ctype.push_back(HEX);
        c2v_tmp.resize(64);
        for (int i = 0; i < c2nv.back(); i++)
          meshFile >> c2v_tmp[i];
        break;

      case 93:
        // Quartic Hexahedron
        c2nv.push_back(125);
        c2nf.push_back(6);
        ctype.push_back(HEX);
        c2v_tmp.resize(125);
        for (int i = 0; i < c2nv.back(); i++)
          meshFile >> c2v_tmp[i];
        break;

      case 94:
        // Quintic Hexahedron
        c2nv.push_back(216);
        c2nf.push_back(6);
        ctype.push_back(HEX);
        c2v_tmp.resize(216);
        for (int i = 0; i < c2nv.back(); i++)
          meshFile >> c2v_tmp[i];
        break;

      default:
        cout << "Gmsh element ID " << k << ", Gmsh Element Type = " << eType
             << endl;
        fatalError("element type not recognized");
        break;
      }

      // Increase the size of c2v (max # of vertices per cell) if needed
      if (c2v.getDim1() < (uint)c2nv[ic]) {
        for (int dim = c2v.getDim1(); dim < c2nv[ic]; dim++) {
          c2v.addCol();
        }
      }

      // Number of nodes in c2v_tmp may vary, so use pointer rather than vector
      c2v.insertRow(c2v_tmp.data(), -1, c2nv[ic]);

      // Shift every value of c2v by -1 (Gmsh is 1-indexed; we need 0-indexed)
      for (int k = 0; k < c2nv[ic]; k++) {
        if (c2v(ic, k) != 0) {
          c2v(ic, k)--;
        }
      }

      ic++;
      getline(meshFile, str); // skip end of line
    } else {
      // Boundary cell; put vertices into bndPts
      int nPtsFace = 0;
      switch (eType) {
      case 1: // Linear edge
        nPtsFace = 2;
        break;

      case 2: // Linear triangle
        nPtsFace = 3;
        break;

      case 3:  // Linear quad
      case 10: // Quadratic (Lagrange) quad
      case 16: // Quadratic (Serendipity) quad
      case 36: // Cubic quad
      case 37: // Quartic quad
      case 38: // Quintic quad
        nPtsFace = 4;
        break;

      case 8: // Quadratic edge
        nPtsFace = 3;
        break;

      case 26: // Cubic Edge
        nPtsFace = 4;
        break;

      case 27: // Quartic Edge
        nPtsFace = 5;
        break;

      case 28: // Quintic Edge
        nPtsFace = 6;
        break;

      case 62: // Order 6
        nPtsFace = 7;
        break;

      case 63: // Order 7
        nPtsFace = 8;
        break;

      case 64: // Order 8
        nPtsFace = 9;
        break;

      case 65: // Order 9
        nPtsFace = 10;
        break;

      case 66: // Order 10
        nPtsFace = 11;
        break;

      default:
        cout << "Gmsh element ID " << k << ", Gmsh Element Type = " << eType
             << endl;
        fatalError("Boundary Element (Face) Type Not Recognized!");
      }

      for (int i = 0; i < nPtsFace; i++) {
        uint32_t iv;
        meshFile >> iv;
        iv--;
        boundPoints[bcid].insert(iv);
      }
      getline(meshFile, str);
    }
  } // End of loop over entities

  int maxNBndPts = 0;
  for (int i = 0; i < nBounds; i++) {
    nBndPts[i] = boundPoints[i].size();
    maxNBndPts = max(maxNBndPts, nBndPts[i]);
  }

  // Copy temp boundPoints data into bndPts matrix
  bndPts.setup(nBounds, maxNBndPts);
  for (int i = 0; i < nBounds; i++) {
    int j = 0;
    for (auto &it : boundPoints[i]) {
      bndPts(i, j) = it;
      j++;
    }
  }

  nEles = c2v.getDim0();

  meshFile.close();
}

bool geo::compareFaces(vector<int> &face1, vector<int> &face2) {
  uint nv = face1.size();
  if (face2.size() != nv)
    return false;

  std::sort(face1.begin(), face1.end());
  std::sort(face2.begin(), face2.end());

  bool found = true;
  for (uint i = 0; i < nv; i++) {
    if (face1[i] != face2[i])
      found = false;
  }

  return found;
}

int geo::compareOrientation(int ic1, int f1, int ic2, int f2) {
  if (nDims == 2)
    return 1;

  int nv = f2nv[c2f(ic1, f1)];

  vector<int> tmpFace1(nv), tmpFace2(nv);

  switch (ctype[ic1]) {
  case HEX:
    // Flux points arranged in 2D grid on each face oriented with each
    // dimension increasing in its +'ve direction ['btm-left' to 'top-right']
    // Node ordering reflects this: CCW from 'bottom-left' node on each face
    switch (f1) {
    case 0:
      // Bottom face  (z = -1)
      tmpFace1[0] = c2v(ic1, 0);
      tmpFace1[1] = c2v(ic1, 1);
      tmpFace1[2] = c2v(ic1, 2);
      tmpFace1[3] = c2v(ic1, 3);
      break;
    case 1:
      // Top face  (z = +1)
      tmpFace1[0] = c2v(ic1, 5);
      tmpFace1[1] = c2v(ic1, 4);
      tmpFace1[2] = c2v(ic1, 7);
      tmpFace1[3] = c2v(ic1, 6);
      break;
    case 2:
      // Left face  (x = -1)
      tmpFace1[0] = c2v(ic1, 0);
      tmpFace1[1] = c2v(ic1, 3);
      tmpFace1[2] = c2v(ic1, 7);
      tmpFace1[3] = c2v(ic1, 4);
      break;
    case 3:
      // Right face  (x = +1)
      tmpFace1[0] = c2v(ic1, 2);
      tmpFace1[1] = c2v(ic1, 1);
      tmpFace1[2] = c2v(ic1, 5);
      tmpFace1[3] = c2v(ic1, 6);
      break;
    case 4:
      // Front face  (y = -1)
      tmpFace1[0] = c2v(ic1, 1);
      tmpFace1[1] = c2v(ic1, 0);
      tmpFace1[2] = c2v(ic1, 4);
      tmpFace1[3] = c2v(ic1, 5);
      break;
    case 5:
      // Back face  (y = +1)
      tmpFace1[0] = c2v(ic1, 3);
      tmpFace1[1] = c2v(ic1, 2);
      tmpFace1[2] = c2v(ic1, 6);
      tmpFace1[3] = c2v(ic1, 7);
      break;
    }
    break;

  default:
    fatalError("Element type not supported.");
    break;
  }

  switch (ctype[ic2]) {
  case HEX:
    switch (f2) {
    case 0:
      // Bottom face  (z = -1)
      tmpFace2[0] = c2v(ic2, 0);
      tmpFace2[1] = c2v(ic2, 1);
      tmpFace2[2] = c2v(ic2, 2);
      tmpFace2[3] = c2v(ic2, 3);
      break;
    case 1:
      // Top face  (z = +1)
      tmpFace2[0] = c2v(ic2, 5);
      tmpFace2[1] = c2v(ic2, 4);
      tmpFace2[2] = c2v(ic2, 7);
      tmpFace2[3] = c2v(ic2, 6);
      break;
    case 2:
      // Left face  (x = -1)
      tmpFace2[0] = c2v(ic2, 0);
      tmpFace2[1] = c2v(ic2, 3);
      tmpFace2[2] = c2v(ic2, 7);
      tmpFace2[3] = c2v(ic2, 4);
      break;
    case 3:
      // Right face  (x = +1)
      tmpFace2[0] = c2v(ic2, 2);
      tmpFace2[1] = c2v(ic2, 1);
      tmpFace2[2] = c2v(ic2, 5);
      tmpFace2[3] = c2v(ic2, 6);
      break;
    case 4:
      // Front face  (y = -1)
      tmpFace2[0] = c2v(ic2, 1);
      tmpFace2[1] = c2v(ic2, 0);
      tmpFace2[2] = c2v(ic2, 4);
      tmpFace2[3] = c2v(ic2, 5);
      break;
    case 5:
      // Back face  (y = +1)
      tmpFace2[0] = c2v(ic2, 3);
      tmpFace2[1] = c2v(ic2, 2);
      tmpFace2[2] = c2v(ic2, 6);
      tmpFace2[3] = c2v(ic2, 7);
      break;
    }
    break;

  default:
    fatalError("Element type not supported.");
    break;
  }

  // Now, compare the two faces to see the relative orientation [rotation]
  if (tmpFace1[0] == tmpFace2[0])
    return 0;
  else if (tmpFace1[1] == tmpFace2[0])
    return 1;
  else if (tmpFace1[2] == tmpFace2[0])
    return 2;
  else if (tmpFace1[3] == tmpFace2[0])
    return 3;
  else
    fatalError("Internal faces improperly matched.");
}