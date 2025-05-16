
#pragma once

#include "ele.hpp"
#include "geo.hpp"
#include "global.hpp"
#include "solver.hpp"

/*! Write solution to file (of type params->plotType) */
void writeData(solver *Solver, input *params);

/*! Write solution data to a CSV file. */
void writeCSV(solver *Solver, input *params);

/*! Write solution data to a Paraview .vtu file. */
void writeParaview(solver *Solver, input *params);

/*! Write out surface data to a Paraview .vtu file. */
void writeSurfaces(solver *Solver, input *params);

/*! Compute the residual and print to both the terminal and history file. */
void writeResidual(solver *Solver, input *params);

/*! Compute and display all error norms */
void writeAllError(solver *Solver, input *params);

/*! Compute the conservation, L1, or L2 solution error (for certain test cases)
 * and print to screen. */
void writeError(solver *Solver, input *params);

/*! Write a Tecplot mesh file compatible with TIOGA's testTioga FORTRAN
 * interface */
void writeMeshTecplot(geo *Geo, input *params);
