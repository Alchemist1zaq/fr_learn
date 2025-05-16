
#pragma once

#include "global.hpp"

#include "input.hpp"
#include "matrix.hpp"

/*! Calculate the inviscid portion of the Euler or Navier-Stokes flux vector at
 * a point */
void inviscidFlux(const double *U, double F[3][5], input *params);

/*! Calculate the viscous portion of the Navier-Stokes flux vector at a point */
void viscousFlux(double *U, matrix<double> &gradU, double Fvis[3][5],
                 input *params);

/*! Calculate the viscous shear-stress tensor at a point */
matrix<double> viscousStressTensor(double *U, matrix<double> &gradU,
                                   input *params);

/*! Calculate the viscous flux for Advection-Diffusion */
void viscousFluxAD(matrix<double> &gradU, double Fvis[3][5], input *params);

/*! Calculate the common inviscid flux at a point using Roe's method */
void roeFlux(double *uL, double *uR, double *norm, double *Fn, input *params);

/*! Calculate the common inviscid flux at a point using the Rusanov
 * scalar-diffusion method */
void rusanovFlux(double *UL, double *UR, double FL[3][5], double FR[3][5],
                 double *norm, double *Fn, double *waveSp, input *params);

/*! Simple central-difference flux (For advection problems) */
void centralFlux(double *uL, double *uR, double *norm, double *Fn,
                 input *params);

/*! Simple central-difference flux (For Navier-Stokes problems) */
void centralFlux(double FL[3][5], double FR[3][5], double *norm, double *Fn,
                 input *params);

/*! Simple upwinded flux (primarily for advection problems) */
void upwindFlux(double *uL, double *uR, double *norm, double *Fn,
                input *params);

/*! Lax-Friedrichs flux (advection-diffusion) */
void laxFriedrichsFlux(double *uL, double *uR, double *norm, double *Fn,
                       input *params);

/*! Calculate the common viscous flux at a point using the LDG penalty method */
void ldgFlux(double *uL, double *uR, matrix<double> &gradU_L,
             matrix<double> &gradU_R, double *Fn, input *params);
