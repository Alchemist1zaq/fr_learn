
#include "solver.hpp"

/* ---- My New Overset Grid Functions ---- */

void solver::insertElement(uint ele_ind) {
  U_spts.add_dim_1(ele_ind, 0.);
  U_fpts.add_dim_1(ele_ind, 0.);

  F_spts.add_dim_2(ele_ind, 0.);
  F_fpts.add_dim_2(ele_ind, 0.);

  V_spts.add_dim_1(ele_ind, 0.);
  V_ppts.add_dim_1(ele_ind, 0.);

  if (params->viscous || params->motion) {
    dU_spts.add_dim_2(ele_ind, 0.);
    dU_fpts.add_dim_2(ele_ind, 0.);
  }

  if (params->viscous) {
    dUc_fpts.add_dim_1(ele_ind, 0.);
  }

  for (auto &mat : dF_spts.data)
    mat.add_dim_1(ele_ind, 0.);

  U0.add_dim_1(ele_ind, 0.);
  U_mpts.add_dim_1(ele_ind, 0.);

  disFn_fpts.add_dim_1(ele_ind, 0.);
  Fn_fpts.add_dim_1(ele_ind, 0.);

  for (auto &divF : divF_spts)
    divF.add_dim_1(ele_ind, 0.);

  /* Multigrid Variables */
  if (params->PMG) {
    sol_spts.add_dim_1(ele_ind, 0.);
    corr_spts.add_dim_1(ele_ind, 0.);
    src_spts.add_dim_1(ele_ind, 0.);
  }

  tempVars_spts.add_dim_1(ele_ind, 0.);
  tempVars_fpts.add_dim_1(ele_ind, 0.);

  pos_spts.add_dim_1(ele_ind, 0.);
  pos_fpts.add_dim_1(ele_ind, 0.);
  pos_ppts.add_dim_1(ele_ind, 0.);

  Jac_spts.add_dim_2(ele_ind, 0.);
  Jac_fpts.add_dim_2(ele_ind, 0.);
  JGinv_spts.add_dim_2(ele_ind, 0.);
  JGinv_fpts.add_dim_2(ele_ind, 0.);
  detJac_spts.add_dim_1(ele_ind, 0.);
  detJac_fpts.add_dim_1(ele_ind, 0.);
  dA_fpts.add_dim_1(ele_ind, 0.);
  norm_fpts.add_dim_1(ele_ind, 0.);

  nodes.add_dim_1(ele_ind, 0.);

  if (params->motion) {
    gridV_spts.add_dim_1(ele_ind, 0.);
    gridV_fpts.add_dim_1(ele_ind, 0.);
    gridV_mpts.add_dim_1(ele_ind, 0.);
    gridV_ppts.add_dim_1(ele_ind, 0.);

    nodesRK.add_dim_1(ele_ind, 0.);
  }

  nEles++;
}

void solver::removeElement(uint ele_ind) {
  U_spts.remove_dim_1(ele_ind);
  U_fpts.remove_dim_1(ele_ind);

  F_spts.remove_dim_2(ele_ind);
  F_fpts.remove_dim_2(ele_ind);

  V_spts.remove_dim_1(ele_ind);
  V_ppts.remove_dim_1(ele_ind);

  if (params->viscous || params->motion) {
    dU_spts.remove_dim_2(ele_ind);
    dU_fpts.remove_dim_2(ele_ind);
  }

  if (params->viscous) {
    dUc_fpts.remove_dim_1(ele_ind);
  }

  for (auto &mat : dF_spts.data)
    mat.remove_dim_1(ele_ind);

  U0.remove_dim_1(ele_ind);
  U_mpts.remove_dim_1(ele_ind);

  disFn_fpts.remove_dim_1(ele_ind);
  Fn_fpts.remove_dim_1(ele_ind);

  for (auto &divF : divF_spts)
    divF.remove_dim_1(ele_ind);

  /* Multigrid Variables */
  if (params->PMG) {
    sol_spts.remove_dim_1(ele_ind);
    corr_spts.remove_dim_1(ele_ind);
    src_spts.remove_dim_1(ele_ind);
  }

  tempVars_spts.remove_dim_1(ele_ind);
  tempVars_fpts.remove_dim_1(ele_ind);

  pos_spts.remove_dim_1(ele_ind);
  pos_fpts.remove_dim_1(ele_ind);
  pos_ppts.remove_dim_1(ele_ind);

  Jac_spts.remove_dim_2(ele_ind);
  Jac_fpts.remove_dim_2(ele_ind);
  JGinv_spts.remove_dim_2(ele_ind);
  JGinv_fpts.remove_dim_2(ele_ind);
  detJac_spts.remove_dim_1(ele_ind);
  detJac_fpts.remove_dim_1(ele_ind);
  dA_fpts.remove_dim_1(ele_ind);
  norm_fpts.remove_dim_1(ele_ind);

  nodes.remove_dim_1(ele_ind);

  if (params->motion) {
    gridV_spts.remove_dim_1(ele_ind);
    gridV_fpts.remove_dim_1(ele_ind);
    gridV_mpts.remove_dim_1(ele_ind);
    gridV_ppts.remove_dim_1(ele_ind);

    nodesRK.remove_dim_1(ele_ind);
  }

  nEles--;
}