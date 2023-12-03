#include "functions.h"

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

void get_exact_soln(int nx, int ny, double *xc, double *yc, double **Tex)
{
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     Tex[i][j] = xc[i]*(1.0-xc[i])*cos(M_PI*yc[j]);
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

