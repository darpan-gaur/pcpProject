#include "utilities.h"
#include "linsolvers.h"

void solve_gssor(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **Tpad, double **Tpnew, int max_iter, double tol, double relax_T)
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

 // now perform iterations
  for(iter=0; iter<max_iter; iter++)
  {
    // update Tpnew
    for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
     {
       ip = i+1; jp = j+1;
       //T_gs = ( b[i][j] + aE[i][j]*Tpad[ip+1][jp] + aW[i][j]*Tpad[ip-1][jp] + 
       //                   aN[i][j]*Tpad[ip][jp+1] + aS[i][j]*Tpad[ip][jp-1] ) / aP[i][j];
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
    printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
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

  // copy from padded array
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
      ip = i+1; jp = j+1;
      T[i][j] = Tpad[ip][jp];
   }

  printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
}

