
#include "multigrid.hpp"

#include <omp.h>
#include <sstream>

#include "cblas.h"

#include "input.hpp"
#include "output.hpp"
#include "solver.hpp"

void multiGrid::setup(int order, input *params, solver &Solver) {
  this->order = order;
  this->params = params;

  pInputs.assign(order, *params);
  pGrids.resize(order);

  if (params->HMG && params->n_h_levels == 0)
    params->HMG = 0;

  /* H-P Multigrid using Refinement Method */
  if (params->HMG) {
    if (params->lowOrder != 0)
      fatalError("H-Multigrid only supported for PMG lowOrder = 0.");

    if (params->nDims == 3)
      fatalError("H-Multigrid only supported for 2D currently.");

    geo coarse_grid;
    coarse_grid.setup(params, true);

    hInputs.assign(params->n_h_levels, *params);
    hGrids.resize(params->n_h_levels);
    hGeos.resize(params->n_h_levels);
    pGeos.resize(params->order);

    for (int H = 0; H < params->n_h_levels; H++) {
      if (params->rank == 0)
        cout << endl << "H-Multigrid: Setting up H = " << H << endl;

      hInputs[H].dataFileName += "_H" + std::to_string((long long)H) + "_";
      hGrids[H] = make_shared<solver>();
      hGeos[H] = make_shared<geo>();

      /* Refine the initial coarse grid to produce the fine grids */
      setup_h_level(coarse_grid, *hGeos[H], params->n_h_levels - H - 1);

      hGrids[H]->setup(&hInputs[H], params->lowOrder, &(*hGeos[H]));
      hGrids[H]->initializeSolution(true);
    }

    /* Create final fine grid and re-setup the given solver */
    if (params->rank == 0)
      cout << endl << "H-Multigrid: Setting up fine grid solver" << endl;
    fine_grid = make_shared<geo>();
    setup_h_level(coarse_grid, *fine_grid, params->n_h_levels);
    Solver.setup(params, order, &(*fine_grid));
    Solver.initializeSolution(false);

    /* Instantiate P-grid solvers using finest mesh */
    for (int P = 0; P < order; P++) {
      if (P < params->lowOrder) {
        // pGrids[P] = NULL;
      } else {
        if (params->rank == 0)
          cout << endl << "P-Multigrid: Setting up P = " << P << endl;

        pInputs[P].dataFileName += "_P" + std::to_string((long long)P) + "_";
        pGrids[P] = make_shared<solver>();
        pGeos[P] = make_shared<geo>();

        setup_h_level(coarse_grid, *pGeos[P], params->n_h_levels);

        pGrids[P]->setup(&pInputs[P], P, &(*pGeos[P]));
        pGrids[P]->initializeSolution(true);
      }
    }
  }

  /* P-Multigrid Alone */
  else {
    /* Instantiate coarse grid solvers */
    for (int P = 0; P < order; P++) {
      if (P < params->lowOrder) {
        // pGrids[P] = NULL;
      } else {
        if (params->rank == 0)
          cout << endl << "P-Multigrid: Setting up P = " << P << endl;

        pInputs[P].dataFileName += "_P" + std::to_string((long long)P) + "_";
        pInputs[P].order = P;
        pGrids[P] = make_shared<solver>();
        pGrids[P]->setup(&pInputs[P], P);
        pGrids[P]->initializeSolution(true);
      }
    }

    if (params->rank == 0)
      cout << endl << "P-Multigrid: Setting up P = " << params->order << endl;
    Solver.setup(params, params->order);
    Solver.initializeSolution(false);
  }
}

void multiGrid::setup_h_level(geo &mesh_c, geo &mesh_f, int refine_level) {
  if (refine_level > 0)
    refineGrid2D(mesh_c, mesh_f, refine_level, mesh_c.nNodesPerCell,
                 params->shapeOrder);
  else
    mesh_f = mesh_c;

  int nSplit = 1 << params->nDims;
  nSplit = std::pow(nSplit, refine_level);

  mesh_f.setup_hmg(params, mesh_c.gridID, mesh_c.gridRank, mesh_c.nProcGrid);
}

void multiGrid::cycle(solver &Solver) {
  /* Update residual on finest grid level and restrict */
  Solver.calcResidual(0);

  restrict_pmg(Solver, *pGrids[order - 1]);

  for (int P = order - 1; P >= (int)params->lowOrder; P--) {
    /* Generate source term */
    compute_source_term(*pGrids[P]);

    /* Copy initial solution to solution storage */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < pGrids[P]->nSpts; spt++)
      for (uint e = 0; e < pGrids[P]->nEles; e++)
        for (uint k = 0; k < params->nFields; k++)
          pGrids[P]->sol_spts(spt, e, k) = pGrids[P]->U_spts(spt, e, k);

    /* Update solution on coarse level */
    for (uint step = 0; step < params->smoothSteps; step++) {
      pGrids[P]->update(true);
    }

    if (P - 1 >= (int)params->lowOrder || params->HMG) {
      /* Update residual and add source */
      pGrids[P]->calcResidual(0);

#pragma omp parallel for collapse(3)
      for (uint spt = 0; spt < pGrids[P]->nSpts; spt++)
        for (uint e = 0; e < pGrids[P]->nEles; e++)
          for (uint k = 0; k < params->nFields; k++)
            pGrids[P]->divF_spts[0](spt, e, k) +=
                pGrids[P]->src_spts(spt, e, k);

      if (P - 1 >= (int)params->lowOrder) {
        /* Restrict to next coarse grid */
        restrict_pmg(*pGrids[P], *pGrids[P - 1]);
      }
    }
  }

  if (params->HMG) {
    /* Initial restriction to first HMG level */
    restrict_hmg(*pGrids[params->lowOrder], *hGrids[0]);

    for (int H = 0; H < params->n_h_levels; H++) {
      /* Generate source term */
      compute_source_term(*hGrids[H]);

      /* Copy initial solution to solution storage */
#pragma omp parallel for collapse(3)
      for (uint spt = 0; spt < hGrids[H]->nSpts; spt++)
        for (uint e = 0; e < hGrids[H]->nEles; e++)
          for (uint k = 0; k < params->nFields; k++)
            hGrids[H]->sol_spts(spt, e, k) = hGrids[H]->U_spts(spt, e, k);

      /* Update solution on coarse level */
      for (uint step = 0; step < params->smoothSteps; step++) {
        hGrids[H]->update(true);
      }

      if (H + 1 < params->n_h_levels) {
        /* Update residual and add source */
        hGrids[H]->calcResidual(0);

#pragma omp parallel for collapse(3)
        for (uint spt = 0; spt < hGrids[H]->nSpts; spt++)
          for (uint e = 0; e < hGrids[H]->nEles; e++)
            for (uint k = 0; k < params->nFields; k++)
              hGrids[H]->divF_spts[0](spt, e, k) +=
                  hGrids[H]->src_spts(spt, e, k);

        /* Restrict to next coarse grid */
        restrict_hmg(*hGrids[H], *hGrids[H + 1]);
      }
    }

    /* --- Upward HMG Cycle --- */
    for (int H = params->n_h_levels - 1; H >= 0; H--) {
      /* Advance again (v-cycle)*/
      for (unsigned int step = 0; step < params->smoothSteps; step++) {
        hGrids[H]->update(true);
      }

      /* Generate error */
#pragma omp parallel for collapse(3)
      for (uint spt = 0; spt < hGrids[H]->nSpts; spt++) {
        for (uint e = 0; e < hGrids[H]->nEles; e++) {
          for (uint k = 0; k < params->nFields; k++) {
            hGrids[H]->corr_spts(spt, e, k) = hGrids[H]->U_spts(spt, e, k);
            hGrids[H]->corr_spts(spt, e, k) -= hGrids[H]->sol_spts(spt, e, k);
          }
        }
      }

      /* Prolong error and add to fine grid solution */
      if (H > 0) {
        prolong_hmg(*hGrids[H], *hGrids[H - 1]);
      } else {
        prolong_hmg(*hGrids[0], *pGrids[params->lowOrder]);
      }
    }
  }

  /* --- Upward PMG Cycle --- */
  for (int P = (int)params->lowOrder; P <= order - 1; P++) {
    /* Advance again (v-cycle)*/
    for (unsigned int step = 0; step < params->smoothSteps; step++) {
      pGrids[P]->update(true);
    }

    /* Generate error */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < pGrids[P]->nSpts; spt++) {
      for (uint e = 0; e < pGrids[P]->nEles; e++) {
        for (uint k = 0; k < params->nFields; k++) {
          pGrids[P]->corr_spts(spt, e, k) = pGrids[P]->U_spts(spt, e, k);
          pGrids[P]->corr_spts(spt, e, k) -= pGrids[P]->sol_spts(spt, e, k);
        }
      }
    }

    /* Prolong error and add to fine grid solution */
    if (P < order - 1) {
      prolong_err(*pGrids[P], *pGrids[P + 1]);
    }
  }

  /* Prolong correction and add to finest grid solution */
  prolong_err(*pGrids[order - 1], Solver);
}

void multiGrid::restrict_pmg(solver &grid_f, solver &grid_c) {
  if (grid_f.order - grid_c.order > 1)
    fatalError("Cannot restrict more than 1 order currently!");

  int m = grid_c.nSpts;
  int k = grid_f.nSpts;
  int n = grid_c.nEles * grid_c.nFields;

  auto &opp_res = grid_f.opers[grid_f.order].opp_restrict(0, 0);

  auto &UF = grid_f.U_spts(0, 0, 0);
  auto &UC = grid_c.U_spts(0, 0, 0);

  auto &dfF = grid_f.divF_spts[0](0, 0, 0);
  auto &dfC = grid_c.divF_spts[0](0, 0, 0);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &opp_res,
              k, &UF, n, 0.0, &UC, n);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &opp_res,
              k, &dfF, n, 0.0, &dfC, n);
}

void multiGrid::prolong_err(solver &grid_c, solver &grid_f) {
  int m = grid_f.nSpts;
  int k = grid_c.nSpts;
  int n = grid_c.nEles * grid_c.nFields;

  auto &opp_pro = grid_c.opers[grid_c.order].opp_prolong(0, 0);
  auto &corr = grid_c.corr_spts(0, 0, 0);
  auto &U = grid_f.U_spts(0, 0, 0);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &opp_pro,
              k, &corr, n, 1.0, &U, n);
}

void multiGrid::compute_source_term(solver &grid) {
  /* Copy restricted fine grid residual to source term */
#pragma omp parallel for collapse(3)
  for (uint spt = 0; spt < grid.nSpts; spt++)
    for (uint e = 0; e < grid.nEles; e++)
      for (uint k = 0; k < params->nFields; k++)
        grid.src_spts(spt, e, k) = grid.divF_spts[0](spt, e, k);

  /* Update residual on current coarse grid */
  grid.calcResidual(0);

  /* Subtract to generate source term */
#pragma omp parallel for collapse(3)
  for (uint spt = 0; spt < grid.nSpts; spt++)
    for (uint e = 0; e < grid.nEles; e++)
      for (uint k = 0; k < params->nFields; k++)
        grid.src_spts(spt, e, k) -= grid.divF_spts[0](spt, e, k);
}

void multiGrid::restrict_hmg(solver &grid_f, solver &grid_c) {
  int nSplit = 1 << params->nDims;

  // what to do...? Take avg (Josh's simple method), or do Galerkin projection?
  // For the moment: take avg value [Only remotely reasonable for P = 0]

  for (uint spt = 0; spt < grid_c.nSpts; spt++) {
    for (uint ec = 0; ec < grid_c.eles.size(); ec++) {

      for (uint k = 0; k < grid_c.nFields; k++) {
        grid_c.U_spts(spt, ec, k) = 0;
        grid_c.divF_spts[0](spt, ec, k) = 0;
      }

      double vol = 0;
      for (uint j = 0; j < nSplit; j++) {
        uint ef = ec * nSplit + j;
        vol += grid_f.detJac_spts(spt, ef);
        for (uint k = 0; k < grid_c.nFields; k++) {
          grid_c.U_spts(spt, ec, k) +=
              grid_f.U_spts(spt, ef, k) * grid_f.detJac_spts(spt, ef);
          grid_c.divF_spts[0](spt, ec, k) +=
              grid_f.divF_spts[0](spt, ef, k) *
              grid_f.detJac_spts(
                  spt,
                  ef); // / nSplit <-- Jameson doesn't average divF, he sums
        }
      }

      for (uint k = 0; k < grid_c.nFields; k++) {
        grid_c.U_spts(spt, ec, k) /= vol;
        grid_c.divF_spts[0](spt, ec, k) /= vol;
      }
    }
  }
}

void multiGrid::prolong_hmg(solver &grid_c, solver &grid_f) {
  int nSplit = 1 << params->nDims;

  // what to do...? Copy value (Josh's simple method), or do Galerkin
  // projection? For the moment: simple method [Only remotely reasonable for P =
  // 0]

  for (uint spt = 0; spt < grid_c.nSpts; spt++) {
    for (uint ef = 0; ef < grid_f.eles.size(); ef++) {
      uint ec = ef / nSplit;
      for (uint k = 0; k < params->nFields; k++) {
        grid_f.U_spts(spt, ef, k) += grid_c.corr_spts(spt, ec, k);
      }
    }
  }
}
