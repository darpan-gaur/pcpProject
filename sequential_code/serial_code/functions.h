#include <math.h>
#include <stdio.h>

void grid(int nx, double xst, double xen, double *xc, double *xf, double *dxc, double *dxf);

void set_initial_guess(int nx, int ny, double *xc, double *yc, double **T);

void calc_diffusivity(int nx, int ny, double *xc, double *yc, double **T, double **kdiff);

void calc_sources(int nx, int ny, double *xc, double *yc, double *dxf, double *dyf, double **T, double **b, double **Sp);

void set_boundary_conditions(int nx, int ny, double *xc, double *yc, double *xf, double *yf, double *dxc, double *dyc, double *dxf, double *dyf, double **T, double *Tleft, double *Tright, double *Ttop, double *qbot);

void get_exact_soln(int nx, int ny, double *xc, double *yc, double **Tex);

void get_fv_coeffs(int nx, int ny, double *xc, double *yc, double *xf, double *yf, double *dxc, double *dxf, double *dyc, double *dyf, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **Sp, double **kdiff, double **T, double *Tleft, double *Tright, double *Ttop, double *qbot);
