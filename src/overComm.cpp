#include "overComm.hpp"

#include "flux.hpp"
#include "global.hpp"

overComm::overComm() {}

void overComm::setup(input *_params, int _nGrids, int _gridID, int _gridRank,
                     int _nprocPerGrid, vector<int> &_gridIdList) {
  params = _params;

  nDims = params->nDims;
  nGrids = _nGrids;
  gridID = _gridID;
  gridRank = _gridRank;
  nprocPerGrid = _nprocPerGrid;
  gridIdList = _gridIdList;

  rank = params->rank;
  nproc = params->nproc;
  nFields = params->nFields;
}

void overComm::setupOverFacePoints(vector<shared_ptr<overFace>> &overFaces,
                                   int nFptsPerFace) {
  // Get all of the fringe points on this grid
  overPts.setup(overFaces.size() * nFptsPerFace, 3);
  overNorm.setup(overFaces.size() * nFptsPerFace, 3);
  uint ipt = 0, jpt = 0;
  for (auto &oface : overFaces) {
    oface->OComm = this;
    oface->fptOffset = ipt;
    auto pts = oface->getPosFpts();
    for (uint i = 0; i < pts.size(); i++) {
      overPts(ipt, 0) = pts[i].x;
      overPts(ipt, 1) = pts[i].y;
      overPts(ipt, 2) = pts[i].z;
      ipt++;
    }
    if (params->oversetMethod == 1) {
      auto norms = oface->getNormFpts();
      for (uint i = 0; i < norms.size(); i++) {
        overNorm(jpt, 0) = norms[i].x;
        overNorm(jpt, 1) = norms[i].y;
        overNorm(jpt, 2) = norms[i].z;
        jpt++;
      }
    }
  }
}

void overComm::setupFringeCellPoints(vector<shared_ptr<ele>> &eles,
                                     const unordered_set<int> &fringeCells,
                                     const vector<int> &eleMap) {
  overPts.setup(0, 0);
  for (auto &ie : fringeCells) {
    int ic = eleMap[ie];
    eles[ic]->sptOffset = overPts.getDim0();
    auto pts = eles[ic]->getPosSpts();
    for (auto &pt : pts)
      overPts.insertRow({pt.x, pt.y, pt.z});
  }
}

void overComm::transferEleData(vector<shared_ptr<ele>> &eles,
                               const unordered_set<int> &fringeCells,
                               const vector<int> &eleMap) {
  int nFields = params->nFields;

  for (auto &ie : fringeCells) {
    int ic = eleMap[ie];
    int Row = eles[ic]->sptOffset;
    for (int spt = 0; spt < eles[ic]->nSpts; spt++) {
      for (int k = 0; k < nFields; k++) {
        eles[ic]->U_spts(spt, k) = U_in(Row + spt, k);
      }
    }
  }
}