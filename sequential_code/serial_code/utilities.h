#include<stdio.h>
#include<math.h>

double get_max_of_array(int nx, int ny, double **arr);
double get_l2err_norm(int nx, int ny, double **arr1, double **arr2);
void output_soln(int nx, int ny, int iter, double *x, double *y, double **T, double **Tex);
