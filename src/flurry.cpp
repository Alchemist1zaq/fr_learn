

#include "flurry.hpp"

#include "funcs.hpp"
#include "multigrid.hpp"

int main(int argc, char *argv[]) {
  input params;
  solver Solver;
  multiGrid pmg;

  int rank = 0;
  int nproc = 1;

  params.rank = rank;
  params.nproc = nproc;

  if (argc < 2) {
    fatalError("No input file specified.");
  }

  /* Read input file & set simulation parameters */
  params.readInputFile(argv[1]);

  if (params.PMG) {
    /* Setup the P-Multigrid class if requested */
    pmg.setup(params.order, &params, Solver);
    std::cout << "PMG" << std::endl;
  } else {
    /* Setup the solver, grid, all elements and faces, and all FR operators for
     * computation */
    Solver.setup(&params, params.order);

    /* Apply the initial condition */
    Solver.initializeSolution();
    std::cout << "NO PMG" << std::endl;
  }

  /* Write initial data file */
  writeData(&Solver, &params);

  double maxTime = params.maxTime;
  int initIter = params.initIter;
  int iterMax = params.iterMax;
  int &iter = params.iter;
  iter = initIter;

  /* --- Calculation Loop --- */
  while (params.iter < iterMax and params.time < maxTime) {
    iter++;

    Solver.update();

    /* If using multigrid, perform correction cycle */
    if (params.PMG)
      pmg.cycle(Solver);

    if ((iter) % params.monitorResFreq == 0 or iter == initIter + 1 or
        params.time >= maxTime)
      writeResidual(&Solver, &params);
    if ((iter) % params.monitorErrFreq == 0 or iter == initIter + 1)
      writeError(&Solver, &params);
    if ((iter) % params.plotFreq == 0 or iter == iterMax or
        params.time >= maxTime)
      writeData(&Solver, &params);
  }

  /* Calculate the integral / L1 / L2 error for the final time */
  writeAllError(&Solver, &params);

  return 0;
}
