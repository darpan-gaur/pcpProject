#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
using namespace std;
using namespace std::chrono;
#define PI 4.0 * atan(1.0)

void grid(int nx, double xst, double xen, double *xc, double *xf, double *dxc, double *dxf)
{
  int i;
  double dxunif;

  // uniform mesh for now; 
  // can use stretching factors to place xc nodes later
  dxunif = (xen-xst)/(double)nx;

  // xc[i] s are inside the CVs
  for(i=0; i<nx; i++)
    xc[i] = ((double)i + 0.5) * dxunif;

  // xf[i] s are at the faces; mid-way betwee xc[i]s
  xf[0] = xst;
  for(i=1; i<nx; i++)
    xf[i] = 0.5 * (xc[i] + xc[i-1]);
  xf[nx] = xen;

  // dxc[i] s are spacing between adjacent xcs; ends are different 
  dxc[0] = xc[0]-xf[0];
  for(i=1; i<nx; i++)
    dxc[i] = xc[i] - xc[i-1];
  dxc[nx] = xf[nx]-xc[nx-1];

  // dxf[i] s are spacing between adjacent xes
  for(i=0; i<nx; i++)
    dxf[i] = xf[i+1] - xf[i];

  //// debug -- print xc
  //printf("--xc--\n");
  //for(i=0; i<nx; i++)
  //  printf("%d %lf\n", i, xc[i]);

  //// debug -- print xf
  //printf("--xf--\n");
  //for(i=0; i<nx+1; i++)
  //  printf("%d %lf\n", i, xf[i]);

  //// debug -- print dxc
  //printf("--dxc--\n");
  //for(i=0; i<nx+1; i++)
  //  printf("%d %lf\n", i, dxc[i]);

  //// debug -- print dxf
  //printf("--dxf--\n");
  //for(i=0; i<nx; i++)
  //  printf("%d %lf\n", i, dxf[i]);
}

void set_initial_guess(int nx, int ny, double *xc, double *yc, double **T)
{
  int i, j;

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      T[i][j] = 0.0;
}

void calc_diffusivity(int nx, int ny, double *xc, double *yc, double **T, double **kdiff)
{
  int i, j;

  // calculate diffusivity (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     kdiff[i][j] = 1.0; // set to uniform for now
   }
}

void calc_sources(int nx, int ny, double *xc, double *yc, double *dxf, double *dyf, double **T, double **b, double **Sp)
{
  int i, j;

  // calculate Su source (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     b[i][j]  = (2.0 + M_PI*M_PI*xc[i]*(1-xc[i])) * cos(M_PI*yc[j]);  // set source function here
     b[i][j] *= dxf[i]*dyf[j];                                  // multiply by volume of CV
   }

  // calculate Sp source (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     Sp[i][j]  = 0.0;           // set source function here
     Sp[i][j] *= dxf[i]*dyf[j]; // multiply by volume of CV
   }
}

void set_boundary_conditions(int nx, int ny, double *xc, double *yc, double *xf, double *yf, double *dxc, double *dyc, double *dxf, double *dyf, double **T, double *Tleft, double *Tright, double *Ttop, double *qbot)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks;

  // left boundary -- Homogeneous Dirichlet
  for(j=0; j<ny-1; j++)
    Tleft[j] = 0.0;

  // right boundary -- Homogeneous Dirichlet
  for(j=0; j<ny-1; j++)
    Tright[j] = 0.0;

  // top boundary -- Homogeneous Dirichlet
  for(i=0; i<nx-1; i++)
    Ttop[i] = 0.0;

  // bottom boundary -- Homogeneous Neumann
  for(i=0; i<nx-1; i++)
    qbot[i] = 0.0;
}


void get_fv_coeffs(int nx, int ny, double *xc, double *yc, double *xf, double *yf, double *dxc, double *dxf, double *dyc, double *dyf, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **Sp, double **kdiff, double **T, double *Tleft, double *Tright, double *Ttop, double *qbot)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks;

  // calculate diffusivity at [xc, yc]; may be dependent on T
  calc_diffusivity(nx, ny, xc, yc, T, kdiff);

  // calculate sources Su, Sp at [xc, yc]; may be dependent on T
  calc_sources(nx, ny, xc, yc, dxf, dyf, T, b, Sp);

  // populate values in BC arrays
  set_boundary_conditions(nx, ny, xc, yc, xf, yf, dxc, dyc, dxf, dyf, 
                          T, Tleft, Tright, Ttop, qbot);

  // start populating the coefficients

  // ------ Step 1 :: interior points ------
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);      aE[i][j] = ke * dyf[j] / dxc[i+1];
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);      aW[i][j] = kw * dyf[j] / dxc[i];

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);      aN[i][j] = kn * dxf[i] / dyc[j+1];
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);      aS[i][j] = ks * dxf[i] / dyc[j];

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // source --  has already been populated
   }

  // ------ Step 2  :: 1-boundary points ------
  // ------ Step 2a :: left boundary ----------
  i = 0;
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);      aE[i][j] = ke * dyf[j] / dxc[i+1];
      kw = 0.0;                                    aW[i][j] = 0.0;

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);      aN[i][j] = kn * dxf[i] / dyc[j+1];
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);      aS[i][j] = ks * dxf[i] / dyc[j];

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Dirichlet BC on the left; 
       // set kdiff and boundary coefficient
       kw = 1.0;        bcW = kw * dyf[j] / dxc[i];
       // now modify b, aP as needed
       b[i][j]  += bcW * Tleft[j];        aP[i][j] += bcW;
   }
  // ------ Step 2a :: left boundary done ---

  // ------ Step 2b :: right boundary ----------
  i = nx-1;
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = 0.0;                                    aE[i][j] = 0.0;
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);      aW[i][j] = kw * dyf[j] / dxc[i];

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);      aN[i][j] = kn * dxf[i] / dyc[j+1];
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);      aS[i][j] = ks * dxf[i] / dyc[j];

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Dirichlet BC on the right; 
       // set kdiff and boundary coefficient
       ke = 1.0;        bcE = ke * dyf[j] / dxc[i+1];
       // now modify b, aP as needed
       b[i][j]  += bcE * Tright[j];        aP[i][j] += bcE;
   }
  // ------ Step 2b :: right boundary done ---

  // ------ Step 2c :: bottom boundary ----------
  j = 0;
  for(i=1; i<nx-1; i++)
   {
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);      aE[i][j] = ke * dyf[j] / dxc[i+1];
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);      aW[i][j] = kw * dyf[j] / dxc[i];

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);      aN[i][j] = kn * dxf[i] / dyc[j+1];
      ks = 0.0;                                    aS[i][j] = 0.0;

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Neumann BC on the bottom; 
       // set kdiff and boundary coefficient
       ks = 1.0;        bcS = ks * dyf[j] / dxc[i];
       // now modify b as needed; aP not modified
       b[i][j]  += qbot[i] * dxf[i];
   }
  // ------ Step 2c :: bottom boundary done ---

  // ------ Step 2d :: top boundary ----------
  j = ny-1;
  for(i=1; i<nx-1; i++)
   {
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);      aE[i][j] = ke * dyf[j] / dxc[i+1];
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);      aW[i][j] = kw * dyf[j] / dxc[i];

      // north-south
      kn = 0.0;                                    aN[i][j] = 0.0;
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);      aS[i][j] = ks * dxf[i] / dyc[j];

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Dirichlet BC on the top; 
       // set kdiff and boundary coefficient
       kn = 1.0;        bcN = kn * dxf[i] / dyc[j+1];
       // now modify b, aP as needed
       b[i][j]  += bcN * Ttop[i];        aP[i][j] += bcN;
   }
  // ------ Step 2d :: top boundary done ---

  // ------ Step 3  :: 2-boundary points ------
  // ------ Step 3a :: top-left boundary ----------
  i = 0; j = ny-1;
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);      aE[i][j] = ke * dyf[j] / dxc[i+1];
      kw = 0.0;                                    aW[i][j] = 0.0;

      // north-south
      kn = 0.0;                                    aN[i][j] = 0.0;
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);      aS[i][j] = ks * dxf[i] / dyc[j];

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Dirichlet BC on the left; 
       kw = 1.0;        bcW = kw * dyf[j] / dxc[i];
       b[i][j]  += bcW * Tleft[j];        aP[i][j] += bcW;

      // Dirichlet BC on the top; 
       kn = 1.0;        bcN = kn * dxf[i] / dyc[j+1];
       b[i][j]  += bcN * Ttop[i];        aP[i][j] += bcN;
  // ------ Step 3a :: top-left boundary done ---

  // ------ Step 3b :: top-right boundary ----------
  i = nx-1; j = ny-1;
      // east-west
      ke = 0.0;                                    aE[i][j] = 0.0;
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);      aW[i][j] = kw * dyf[j] / dxc[i];

      // north-south
      kn = 0.0;                                    aN[i][j] = 0.0;
      ks = 0.5*(kdiff[i][j] + kdiff[i][j-1]);      aS[i][j] = ks * dxf[i] / dyc[j];

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Dirichlet BC on the right; 
       ke = 1.0;        bcE = ke * dyf[j] / dxc[i+1];
       b[i][j]  += bcE * Tright[j];        aP[i][j] += bcE;

      // Dirichlet BC on the top; 
       kn = 1.0;        bcN = kn * dxf[j] / dyc[j+1];
       b[i][j]  += bcN * Ttop[j];        aP[i][j] += bcN;
  // ------ Step 3b :: top-right boundary done ---

  // ------ Step 3c :: bottom-left boundary ----------
  i = 0; j = 0;
      // east-west
      ke = 0.5*(kdiff[i][j] + kdiff[i+1][j]);      aE[i][j] = ke * dyf[j] / dxc[i+1];
      kw = 0.0;                                    aW[i][j] = 0.0;

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);      aN[i][j] = kn * dxf[i] / dyc[j+1];
      ks = 0.0;                                    aS[i][j] = 0.0;

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Neumann BC on the bottom; 
       b[i][j]  += qbot[i] * dxf[i];

      // Dirichlet BC on the left; 
       kw = 1.0;        bcW = kw * dyf[j] / dxc[i];
       b[i][j]  += bcW * Tleft[j];        aP[i][j] += bcW;
  // ------ Step 3c :: bottom-left boundary done ---

  // ------ Step 3d :: bottom-right boundary ----------
  i = nx-1;  j = 0;
      // east-west
      ke = 0.0;                                    aE[i][j] = 0.0;
      kw = 0.5*(kdiff[i][j] + kdiff[i-1][j]);      aW[i][j] = kw * dyf[j] / dxc[i];

      // north-south
      kn = 0.5*(kdiff[i][j] + kdiff[i][j+1]);      aN[i][j] = kn * dxf[i] / dyc[j+1];
      ks = 0.0;                                    aS[i][j] = 0.0;

      // present
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j];

      // Neumann BC on the bottom; 
       b[i][j]  += qbot[i] * dxf[i];

      // Dirichlet BC on the right; 
       ke = 1.0;        bcE = ke * dyf[j] / dxc[i+1];
       b[i][j]  += bcE * Tright[j];        aP[i][j] += bcE;
  // ------ Step 3d :: bottom-right boundary done ---

  //// debug
  //for(i=0; i<nx; i++)
  // for(j=0; j<ny; j++)
  // {
  //   printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aP[i][j], aE[i][j], aW[i][j], aN[i][j], aS[i][j], b[i][j]);
  // }
 
}

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

long long solve_gssor(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **Tpad, double **Tpnew, int max_iter, double tol, double relax_T)
{

  int i, j, ip, jp, iter;
  double l2err, arrmax1, arrmax2, err_ref, T_gs, rel_err;

  // copy to padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      Tpad[ip][jp] = T[i][j];
   }

    auto startTime = chrono::high_resolution_clock::now();
 // now perform iterations
  for(iter=0; iter<max_iter; iter++)
  {
    // update Tpnew
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       T_gs = ( b[i][j] + aE[i][j]*Tpnew[ip+1][jp] + aW[i][j]*Tpnew[ip-1][jp] + 
                          aN[i][j]*Tpnew[ip][jp+1] + aS[i][j]*Tpnew[ip][jp-1] ) / aP[i][j];
       Tpnew[ip][jp] = (1.0-relax_T)*Tpad[ip][jp] + relax_T*T_gs;
     }

    // check for convergence
    l2err = get_l2err_norm(nx+2, ny+2, Tpad, Tpnew);
    arrmax1 = get_max_of_array(nx+2, ny+2, Tpad );
    arrmax2 = get_max_of_array(nx+2, ny+2, Tpnew);
    err_ref = fmax(arrmax1, arrmax2);  err_ref = fmax(err_ref, 1.0e-6);
    rel_err = l2err/err_ref;
    if(rel_err < tol)
      break;

    // prepare for next iteration
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       Tpad[ip][jp] = Tpnew[ip][jp];
     }

  }

  auto endTime = chrono::high_resolution_clock::now();
  long long timeTaken = duration_cast<chrono::milliseconds>(endTime - startTime).count();

  // copy from padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      T[i][j] = Tpad[ip][jp];
   }

   printf("\nNumber of iterations: %d\n", min(iter+1, max_iter));
   printf("Final error: %9.5e\n", rel_err);
//   printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
    return timeTaken;
}

void get_exact_soln(int nx, int ny, double *xc, double *yc, double **Tex)
{
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     Tex[i][j] = xc[i]*(1.0-xc[i])*cos(M_PI*yc[j]);
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

int main()
{

  int nx, ny;
  double *xf, *xc, *dxf, *dxc, *yf, *yc, *dyc, *dyf;
  double ** aP, **aE, **aW, **aN, **aS, **b, **Sp;
  double **T, **Tex, **kdiff, **wrk1, **wrk2;
  double xst, xen, yst, yen;
  double *Tleft, *Tright, *Ttop, *qbot;
  int i, j, max_iter;
  double relax_T, tol, l2err;
  FILE* fp;  

  // read inputs
  fp = fopen("input_seq.in", "r");
  fscanf(fp, "%d %d\n", &nx, &ny);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &yst, &yen);
  fclose(fp);

  printf("Inputs are: %d %d %lf %lf %lf %lf\n", nx, ny, xst, xen, yst, yen);

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  xc  = (double *)malloc( nx    * sizeof(double));   // CV centers (interiors)
  xf  = (double *)malloc((nx+1) * sizeof(double));   // CV faces
  dxc = (double *)malloc((nx+1) * sizeof(double));   // spacing betw centers (interiors)
  dxf = (double *)malloc( nx    * sizeof(double));   // spacing betw faces

  yc  = (double *)malloc( ny    * sizeof(double));   // CV centers (interiors)
  yf  = (double *)malloc((ny+1) * sizeof(double));   // CV faces
  dyc = (double *)malloc((ny+1) * sizeof(double));   // spacing betw centers (interiors)
  dyf = (double *)malloc( ny    * sizeof(double));   // spacing betw faces

  Tleft  = (double *)malloc( ny * sizeof(double));   // left   Dirichlet condition
  Tright = (double *)malloc( ny * sizeof(double));   // right  Dirichlet condition
  Ttop   = (double *)malloc( nx * sizeof(double));   // top    Dirichlet condition
  qbot   = (double *)malloc( nx * sizeof(double));   // bottom Neumann   condition
  printf("   >> Done allocating 1D arrays -- \n");


  // allocate 2D arrays dynamically
  // -- for T --
  T = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    T[i] = (double *)malloc(ny*sizeof(double));

  // -- for Tex --
  Tex = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    Tex[i] = (double *)malloc(ny*sizeof(double));

  // -- for aP --
  aP = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aP[i] = (double *)malloc(ny*sizeof(double));

  // -- for aE --
  aE = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aE[i] = (double *)malloc(ny*sizeof(double));

  // -- for aW --
  aW = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aW[i] = (double *)malloc(ny*sizeof(double));

  // -- for aN --
  aN = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aN[i] = (double *)malloc(ny*sizeof(double));

  // -- for aS --
  aS = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aS[i] = (double *)malloc(ny*sizeof(double));

  // -- for b --
  b = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    b[i] = (double *)malloc(ny*sizeof(double));

  // -- for Sp --
  Sp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    Sp[i] = (double *)malloc(ny*sizeof(double));

  // -- for kdiff --
  kdiff = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    kdiff[i] = (double *)malloc(ny*sizeof(double));

  // -- for work arrays wrk1, wrk2 --
  wrk1 = (double **)malloc((nx+2)*sizeof(double *));
  for(i=0; i<nx+2; i++)
    wrk1[i] = (double *)malloc((ny+2)*sizeof(double));

  wrk2 = (double **)malloc((nx+2)*sizeof(double *));
  for(i=0; i<nx+2; i++)
    wrk2[i] = (double *)malloc((ny+2)*sizeof(double));


  printf(" > Done allocating memory -------- \n");
  
  // initialize the grid
  grid(nx, xst, xen, xc, xf, dxc, dxf);  // -- along x --
  grid(ny, yst, yen, yc, yf, dyc, dyf);  // -- along y --
  printf("\n > Done setting up grid ---------- \n");

  set_initial_guess(nx,ny,xc,yc,T);  // initial condition
  printf("\n > Done setting up initial guess -- \n");

  get_fv_coeffs(nx, ny, xc, yc, xf, yf, dxc, dxf, dyc, dyf,     // grid vars
                aP, aE, aW, aN, aS, b, Sp, kdiff, T,            // coefficients
                Tleft, Tright, qbot, Ttop);                     // BC vars
  printf("\n > Done calculating fv coeffs ----- \n");

  printf("\n > Solving for T ------------- \n");
  max_iter = 100000; tol = 1.0e-10; relax_T = 1.0;
  long long timeTaken = solve_gssor(nx, ny, aP, aE, aW, aN, aS, b, T, wrk1, wrk2, max_iter, tol, relax_T);

    printf("Time taken: %lld ms\n\n", timeTaken);
//   printf(" > Done solving for T ------------- \n");

  get_exact_soln(nx, ny, xc, yc, Tex);
  output_soln(nx, ny, 0, xc, yc, T, Tex);

  l2err = get_l2err_norm(nx, ny, T, Tex);
//   printf("%d %d %9.5e", nx, ny, l2err);

  // free memory
   // ----1D arrays ---
   free(yf);  free(yc);    free(dyf);    free(dyc);
   free(xf);  free(xc);    free(dxf);    free(dxc);
   free(Tleft); free(Tright); free(Ttop);  free(qbot);
   // --- Done 1D arrays ---

   // ----2D arrays ---
   for(i=0; i<nx; i++)     free(T[i]);      free(T);
   for(i=0; i<nx; i++)     free(Tex[i]);    free(Tex);
   for(i=0; i<nx; i++)     free(Sp[i]);     free(Sp);
   for(i=0; i<nx; i++)     free(kdiff[i]);  free(kdiff);
   for(i=0; i<nx; i++)     free(b[i]);      free(b);
   for(i=0; i<nx; i++)     free(aP[i]);     free(aP);
   for(i=0; i<nx; i++)     free(aE[i]);     free(aE);
   for(i=0; i<nx; i++)     free(aW[i]);     free(aW);
   for(i=0; i<nx; i++)     free(aN[i]);     free(aN);
   for(i=0; i<nx; i++)     free(aS[i]);     free(aS);
   for(i=0; i<nx+2; i++)   free(wrk1[i]);   free(wrk1);
   for(i=0; i<nx+2; i++)   free(wrk2[i]);   free(wrk2);
   // --- Done 2D arrays ---
//   printf("\n > Done freeing up memory --------- \n");
  return 0;
}