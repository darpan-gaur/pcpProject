#include "utilities.h"

double get_max_of_array(int nx, int ny, double **arr)
{
  int i, j;
  double arrmax, val;

  arrmax = arr[0][0];
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     val = fabs(arr[i][j]);
     if(arrmax < val)
       arrmax = val;
   }
  return arrmax;
}

double get_l2err_norm(int nx, int ny, double **arr1, double **arr2)
{
  double l2err = 0.0, val;
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     val = arr1[i][j] - arr2[i][j];
     l2err += val*val;
   }
   //printf("l2err = %lf\n", l2err);
  l2err = l2err/((double) (nx*ny));
  l2err = sqrt(l2err);

  return l2err;
}

void output_soln(int nx, int ny, int iter, double *x, double *y, double **T, double **Tex)
{
  int i, j;
  FILE* fp;
  char fname[100];

  sprintf(fname, "T_xy_%03d_%03d_%04d.dat", nx, ny, iter);

  fp = fopen(fname, "w");
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     fprintf(fp, "%lf %lf %lf %lf\n", x[i], y[j], T[i][j], Tex[i][j]);
  fclose(fp);

  printf("Done writing solution for stamp = %d to file %s\n", iter, fname);
}

