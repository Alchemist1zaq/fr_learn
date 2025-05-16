
#pragma once

#include <memory>
#include <vector>

#include "input.hpp"
#include "matrix.hpp"
#include "solver.hpp"

class multiGrid {
private:
  input *params = NULL;
  vector<input> pInputs, hInputs;
  int order;
  vector<shared_ptr<solver>> pGrids, hGrids;
  vector<shared_ptr<geo>> pGeos, hGeos;
  shared_ptr<geo> fine_grid;

  //! For nested HMG method
  vector<vector<int>> parent_cells;
  vector<matrix<int>> child_cells;

  void restrict_pmg(solver &grid_fine, solver &grid_coarse);
  void prolong_err(solver &grid_c, solver &grid_f);
  void compute_source_term(solver &grid);

  void restrict_hmg(solver &grid_f, solver &grid_c);
  void prolong_hmg(solver &grid_c, solver &grid_f);
  void setup_h_level(geo &mesh_c, geo &mesh_f, int refine_level);

public:
  void setup(int order, input *params, solver &Solver);
  void cycle(solver &Solver);
};
