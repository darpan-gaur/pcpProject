#include<iostream>
#include<cstring>
#include<thread>
#include<vector>
#include<string>
// #include<barrier>
#include<chrono>
#include<pthread.h>
#include<stdlib.h>
#include<math.h>
#define PI 4.0 * atan(1.0)
using namespace std;
using namespace std::chrono;

// global variables
int nx, ny, num_diagonals, iteration_batch_size;
int num_threads;
double *x, *dx, *y, *dy;
double **aP, **aE, **aW, **aN, **aS, **b, **Sp;
double **T, **Tex, **kdiff, **Tpnew;
double xst, xen, yst, yen;
double **bcleft, **bcrght, **bctop, **bcbot;
int max_iter;
double relax_T, tol, rel_err = 1000.0;
int *thread_ids;
double *thread_errs, *thread_tpad_max, *thread_tpnew_max;
int num_iters_taken;

pthread_barrier_t bar_calc;
pthread_barrier_t bar_upd;

void grid(int nx, double xst, double xen, double *x, double *dx)
{
    int i;
    double dxunif;

    // uniform mesh for now;
    // can use stretching factors to place x nodes later
    dxunif = (xen - xst) / (double)(nx - 1);

    // populate x[i] s
    for (i = 0; i < nx; i++)
        x[i] = xst + ((double)i) * dxunif;

    // dx[i] s are spacing between adjacent xs
    for (i = 0; i < nx - 1; i++)
        dx[i] = dxunif;

    //// debug -- print x
    //   printf("--x--\n");
    //   for(i=0; i<nx; i++)
    //     printf("%d %lf\n", i, x[i]);

    //// debug -- print dx
    //   printf("--dx--\n");
    //   for(i=0; i<nx-1; i++)
    //     printf("%d %lf\n", i, dx[i]);
}

void set_initial_guess(int nx, int ny, double *x, double *y, double **T)
{
    int i, j;

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            T[i][j] = 1.0;
}

void calc_diffusivity(int nx, int ny, double *x, double *y, double **T, double **kdiff)
{
    int i, j;
    double A = 1.0, B = 0.0;

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            kdiff[i][j] = 1;
        }
}

void calc_sources(int nx, int ny, double *x, double *y, double *dx, double *dy, double **T, double **b)
{
    int i, j;

    // calculate source (may be dependent on T)
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            double xval = (double)i * dx[0], yval = (double)j * dy[0];
            b[i][j] = -1.0 * cos(PI * yval) * (2 + PI * PI * xval * (1 - xval)); // set source function here
            b[i][j] = -b[i][j] * dx[0] * dx[0];                                  // only works for uniform mesh
        }
}

void set_boundary_conditions(int nx, int ny, double *x, double *y, double *dx, double *dy, double **T, double **bcleft, double **bcrght, double **bctop, double **bcbot)
{

    int i, j;
    double bcE, bcW, bcN, bcS, kw, ke, kn, ks;

    // left boundary -- Homogeneous Dirichlet
    // q(y) * T + p(y) * dTdx = r(y)
    for (j = 0; j < ny; j++)
    {
        bcleft[j][0] = 0.0; // p(y)
        bcleft[j][1] = 1.0; // q(y)
        bcleft[j][2] = 0.0; // r(y)
    }

    // right boundary -- Homogeneous Dirichlet
    // q(y) * T + p(y) * dTdx = r(y)
    for (j = 0; j < ny; j++)
    {
        bcrght[j][0] = 0.0; // p(y)
        bcrght[j][1] = 1.0; // q(y)
        bcrght[j][2] = 0.0; // r(y)
    }

    // bottom boundary -- Homogeneous  Neumann
    // q(x) * T + p(x) * dTdy = r(x)
    for (i = 0; i < nx; i++)
    {
        bcbot[i][0] = 1.0; // p(x)
        bcbot[i][1] = 0.0; // q(x)
        bcbot[i][2] = 0.0; // r(x)
    }

    // top boundary -- Homogeneous Dirichlet
    // q(x) * T + p(x) * dTdy = r(x)
    for (i = 0; i < nx; i++)
    {
        bctop[i][0] = 0.0; // p(x)
        bctop[i][1] = 1.0; // q(x)
        bctop[i][2] = 0.0; // r(x)
    }
}

void get_coeffs(int nx, int ny, double *x, double *y, double *dx, double *dy, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **T, double **kdiff, double **bcleft, double **bcrght, double **bctop, double **bcbot)
{

    int i, j;
    double pj, qj, rj, pi, qi, ri;
    double hxhy_sq;

    // calculate diffusivity at [x, y]; may be dependent on T
    calc_diffusivity(nx, ny, x, y, T, kdiff);

    // calculate sources Su, Sp at [x, y]; may be dependent on T
    calc_sources(nx, ny, x, y, dx, dy, T, b);

    // populate values in BC arrays
    set_boundary_conditions(nx, ny, x, y, dx, dy, T, bcleft, bcrght, bctop, bcbot);
    //   for (i=0;i<nx;i++) {
    //     printf("Bottom pi[%d] = %lf\n",i, bcbot[i][0]);
    //   }

    // start populating the coefficients
    hxhy_sq = dx[0] * dx[0] / (dy[0] * dy[0]);

    // ------ Step 1 :: interior points ------
    for (i = 1; i < nx - 1; i++)
        for (j = 1; j < ny - 1; j++)
        {
            aP[i][j] = 2.0 * (1.0 + hxhy_sq);
            aE[i][j] = -1.0;
            aW[i][j] = -1.0;
            aN[i][j] = -1.0;
            aS[i][j] = -1.0;
        }

    // ------ Step 2 :: left boundary ----------
    i = 0;
    for (j = 0; j < ny; j++)
    {
        pj = bcleft[j][0];
        qj = bcleft[j][1];
        rj = bcleft[j][2];

        aP[i][j] = 2.0 * (1.0 + hxhy_sq) * pj - 2.0 * dx[0] * qj;
        aE[i][j] = -2.0 * pj;
        aW[i][j] = 0.0;
        aN[i][j] = -hxhy_sq * pj;
        aS[i][j] = -hxhy_sq * pj;
        b[i][j] = b[i][j] * pj - 2.0 * dx[0] * rj;
    }
    // ------ Step 2 :: left boundary done ---

    // ------ Step 3 :: right boundary ----------
    i = nx - 1;
    for (j = 0; j < ny; j++)
    {
        pj = bcrght[j][0];
        qj = bcrght[j][1];
        rj = bcrght[j][2];

        aP[i][j] = 2.0 * (1.0 + hxhy_sq) * pj + 2.0 * dx[0] * qj;
        aE[i][j] = 0.0;
        aW[i][j] = -2.0 * pj;
        aN[i][j] = -hxhy_sq * pj;
        aS[i][j] = -hxhy_sq * pj;
        b[i][j] = b[i][j] * pj + 2.0 * dx[0] * rj;
    }
    // ------ Step 3 :: right boundary done ---

    // ------ Step 4 :: bottom boundary ----------
    j = 0;
    for (i = 1; i < nx - 1; i++)
    {
        pi = bcbot[i][0];
        qi = bcbot[i][1];
        ri = bcbot[i][2];
        //  printf("pi = %f, qi = %f, ri = %f\n", pi, qi, ri);
        aP[i][j] = 2.0 * (1.0 + hxhy_sq) * pi - 2.0 * hxhy_sq * dy[0] * qi;
        //  printf("aP[%d][%d] = %f\n", i, j, aP[i][j]);
        aW[i][j] = -1.0 * pi;
        aE[i][j] = -1.0 * pi;
        aN[i][j] = -2.0 * hxhy_sq * pi;
        aS[i][j] = 0.0;
        b[i][j] = b[i][j] * pi - 2.0 * hxhy_sq * dy[0] * ri;
    }
    // ------ Step 4 :: bottom boundary done ---

    // ------ Step 5 :: top boundary ----------
    j = ny - 1;
    for (i = 1; i < nx - 1; i++)
    {
        pi = bctop[i][0];
        qi = bctop[i][1];
        ri = bctop[i][2];

        aP[i][j] = 2.0 * (1.0 + hxhy_sq) * pi + 2.0 * hxhy_sq * dy[0] * qi;
        aW[i][j] = -1.0 * pi;
        aE[i][j] = -1.0 * pi;
        aN[i][j] = 0.0;
        aS[i][j] = -2.0 * hxhy_sq * pi;
        b[i][j] = b[i][j] * pi + 2.0 * hxhy_sq * dy[0] * ri;
    }
    // ------ Step 5 :: top boundary done ---

    // debug
    //   for(i=0; i<nx; i++)
    //    for(j=0; j<ny; j++)
    //    {
    //      printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aP[i][j], aE[i][j], aW[i][j], aN[i][j], aS[i][j], b[i][j]);
    //    }
}

double get_max_of_array(int nx, int ny, double **arr)
{
    int i, j;
    double arrmax, val;

    arrmax = arr[0][0];
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            val = fabs(arr[i][j]);
            if (arrmax < val)
                arrmax = val;
        }
    return arrmax;
}

double get_max_1d_array(int n, double* arr) {
    int i;
    double arrmax, val;

    arrmax = arr[0];
    for (i = 0; i < n; i++)
    {
        val = fabs(arr[i]);
        if (arrmax < val)
            arrmax = val;
    }
    return arrmax;
}

double get_sum_1d_array(int n, double* arr) {
    double su = 0.0;
    for (int i=0; i<n; i++) {
        su += arr[i];
    }
    return su;
}

double get_l2err_norm(int nx, int ny, double **arr1, double **arr2)
{
    double l2err = 0.0, val;
    int i, j;

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            val = arr1[i][j] - arr2[i][j];
            l2err += val * val;
        }
    // printf("l2err = %lf\n", l2err);
    l2err = l2err / ((double)(nx * ny));
    l2err = sqrt(l2err);

    return l2err;
}

// Runner function of threads which calculate new values of Tpnew
void *calculate_new_values(void* thread_id) {
    int *thrIdPtr = (int *)thread_id;
    int thrId = *thrIdPtr;

    int iter, diagonal_number, ip , jp, i, j, batch_iter_count = 0;
    double T_gs, T_prev, thr_err, diff, tpad_max, tpnew_max;
    for (iter = 0; iter < max_iter; iter++) {
        if (batch_iter_count == iteration_batch_size) {
            // Keep calculating error and max values when batch_iter_count == iteration_batch_size
            thr_err = 0.0;
            tpad_max = 0.0;
            tpnew_max = 0.0;
        }
        // Calculate new values
        diagonal_number = thrId;
        while (diagonal_number <= num_diagonals) {
            // Starting points of the diagonal
            ip = min(diagonal_number, nx);
            jp = max(diagonal_number - nx +1, 1);

            while (ip >= 1 and jp <= ny) {

                i = ip-1; j = jp-1;

                T_prev = Tpnew[ip][jp];
                T_gs = (b[i][j] - aE[i][j] * Tpnew[ip + 1][jp] - aW[i][j] * Tpnew[ip - 1][jp] -
                        aN[i][j] * Tpnew[ip][jp + 1] - aS[i][j] * Tpnew[ip][jp - 1]) /
                       aP[i][j];
                Tpnew[ip][jp] = (1.0 - relax_T) * T_prev + relax_T * T_gs;

                if (batch_iter_count == iteration_batch_size) {
                    // Calculate error
                    diff = Tpnew[ip][jp] - T_prev;
                    thr_err += diff * diff;

                    // Update max values
                    tpad_max = max(tpad_max, fabs(T_prev));
                    tpnew_max = max(tpnew_max, fabs(Tpnew[ip][jp]));
                }

                ip--;
                jp++;
            }
            
            diagonal_number += num_threads;
        }

        if (batch_iter_count == iteration_batch_size) {
            // Update thread error and max values
            thread_errs[thrId-1] = thr_err;
            thread_tpad_max[thrId-1] = tpad_max;
            thread_tpnew_max[thrId-1] = tpnew_max;

            pthread_barrier_wait(&bar_calc);
            pthread_barrier_wait(&bar_upd);

            if (rel_err < tol) break;
        }

        if (batch_iter_count == iteration_batch_size) batch_iter_count = 0;
        else batch_iter_count++;
    }
    num_iters_taken = min(iter+1, max_iter);
    pthread_exit(0);
}

void *update_values(void* thread_id) {

    // cout << "Update thread started" << endl;

    int iter, i, j, ip, jp;
    double arrmax1, arrmax2, err_ref, l2err;
    for (iter = 0; iter < max_iter; iter++) {

        pthread_barrier_wait(&bar_calc);

        l2err = get_sum_1d_array(num_threads, thread_errs);
        l2err = l2err / ((double)(nx * ny));
        l2err = sqrt(l2err);

        arrmax1 = get_max_1d_array(num_threads, thread_tpad_max);
        arrmax2 = get_max_1d_array(num_threads, thread_tpnew_max);
        err_ref = fmax(arrmax1, arrmax2);
        err_ref = fmax(err_ref, 1.0e-6);
        rel_err = l2err / err_ref;
        if (rel_err < tol) {
            pthread_barrier_wait(&bar_upd);
            break;
        }

        pthread_barrier_wait(&bar_upd);
    }
    pthread_exit(0);
}

long long solve_gssor_wave()
{

    int i, j, ip, jp;
    double l2err, arrmax1, arrmax2, err_ref, T_gs;

    // copy to padded array
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++)
        {
            ip = i + 1;
            jp = j + 1;
            Tpnew[ip][jp] = T[i][j];
        }
    }

    auto startTime = chrono::high_resolution_clock::now();

    pthread_t calculate_new_values_threads[num_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);

    for (i=0; i < num_threads; i++) {
        pthread_create(&calculate_new_values_threads[i], &attr, calculate_new_values, thread_ids + i);
    }

    pthread_t update_thread;
    int updatedId = 0;
    pthread_create(&update_thread, &attr, update_values, &updatedId);

    // Join all threads
    for (i=0; i<num_threads; i++) {
        pthread_join(calculate_new_values_threads[i], NULL);
    }

    // update_thread.join();
    pthread_join(update_thread, NULL);

    auto endTime = chrono::high_resolution_clock::now();
    long long timeTaken = duration_cast<chrono::milliseconds>(endTime - startTime).count();

    // copy from padded array
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            ip = i + 1;
            jp = j + 1;
            T[i][j] = Tpnew[ip][jp];
        }

    return timeTaken;
}

void get_exact_soln(int nx, int ny, double *x, double *y, double **Tex)
{
    int i, j;

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            Tex[i][j] = x[i] * (1 - x[i]) * cos(PI * y[j]);
}

void output_soln(int nx, int ny, int iter, double *x, double *y, double **T, double **Tex)
{
    int i, j;
    FILE *fp;
    char fname[100];

    sprintf(fname, "T_xy_par_%03d_%03d_%04d.dat", nx, ny, iter);

    fp = fopen(fname, "w");
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            fprintf(fp, "%lf %lf %lf %lf\n", x[i], y[j], T[i][j], Tex[i][j]);
    fclose(fp);

    printf("Done writing solution for stamp = %d to file %s\n", iter, fname);
}

int main()
{
    int i, j;
    FILE *fp;

    // read inputs
    fp = fopen("input_batch.in", "r");
    fscanf(fp, "%d %d\n", &nx, &ny);
    fscanf(fp, "%lf %lf\n", &xst, &xen);
    fscanf(fp, "%lf %lf\n", &yst, &yen);
    fscanf(fp, "%d %d\n", &num_threads, &iteration_batch_size);
    fclose(fp);

    printf("Inputs are: %d %d %lf %lf %lf %lf %d %d\n", nx, ny, xst, xen, yst, yen, num_threads, iteration_batch_size);

    num_diagonals = nx + ny - 1;

    pthread_barrier_init(&bar_calc, NULL, num_threads+1);
    pthread_barrier_init(&bar_upd, NULL, num_threads+1);

    // allocate memory
    printf("\n > Allocating Memory -- \n");
    x = (double *)malloc(nx * sizeof(double));        // grid points
    dx = (double *)malloc((nx - 1) * sizeof(double)); // spacing betw grid points

    y = (double *)malloc(ny * sizeof(double));        // grid points
    dy = (double *)malloc((ny - 1) * sizeof(double)); // spacing betw grid points

    printf("   >> Done allocating 1D arrays -- \n");

    // allocate 2D arrays dynamically
    // -- for T --
    T = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        T[i] = (double *)malloc(ny * sizeof(double));

    // -- for Tex --
    Tex = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        Tex[i] = (double *)malloc(ny * sizeof(double));

    // -- for aP --
    aP = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        aP[i] = (double *)malloc(ny * sizeof(double));

    // -- for aE --
    aE = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        aE[i] = (double *)malloc(ny * sizeof(double));

    // -- for aW --
    aW = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        aW[i] = (double *)malloc(ny * sizeof(double));

    // -- for aN --
    aN = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        aN[i] = (double *)malloc(ny * sizeof(double));

    // -- for aS --
    aS = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        aS[i] = (double *)malloc(ny * sizeof(double));

    // -- for b --
    b = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        b[i] = (double *)malloc(ny * sizeof(double));

    // left boundary condition
    bcleft = (double **)malloc(ny * sizeof(double *));
    for (i = 0; i < ny; i++)
        bcleft[i] = (double *)malloc(3 * sizeof(double));

    // right boundary condition
    bcrght = (double **)malloc(ny * sizeof(double *));
    for (i = 0; i < ny; i++)
        bcrght[i] = (double *)malloc(3 * sizeof(double));

    // bottom boundary condition
    bctop = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        bctop[i] = (double *)malloc(3 * sizeof(double));

    // top boundary condition
    bcbot = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        bcbot[i] = (double *)malloc(3 * sizeof(double));

    // -- for kdiff --
    kdiff = (double **)malloc(nx * sizeof(double *));
    for (i = 0; i < nx; i++)
        kdiff[i] = (double *)malloc(ny * sizeof(double));

    Tpnew = (double **)malloc((nx + 2) * sizeof(double *));
    for (i = 0; i < nx + 2; i++)
        Tpnew[i] = (double *)malloc((ny + 2) * sizeof(double));

    thread_ids = (int *)malloc(num_threads * sizeof(int));
    thread_errs = (double *)malloc(num_threads * sizeof(double));
    thread_tpad_max = (double *)malloc(num_threads * sizeof(double));
    thread_tpnew_max = (double *)malloc(num_threads * sizeof(double));
    for (i = 0; i < num_threads; i++) {
        thread_ids[i] = i+1;
    }

    printf("   >> Done allocating 2D arrays -- \n");
    printf(" > Done allocating memory -------- \n");

    // initialize the grid
    grid(nx, xst, xen, x, dx); // -- along x --
    grid(ny, yst, yen, y, dy); // -- along y --
    printf("\n > Done setting up grid ---------- \n");

    set_initial_guess(nx, ny, x, y, T); // initial condition
    printf("\n > Done setting up initial guess -- \n");

    // ---
    get_coeffs(nx, ny, x, y, dx, dy,            // grid vars
               aP, aE, aW, aN, aS, b, T, kdiff, // coefficients
               bcleft, bcrght, bcbot, bctop);   // BC vars
    printf("\n > Done calculating coeffs ----- \n");

    printf("\n > Solving for T ------------- \n");
    max_iter = 100000;
    tol = 1.0e-10;
    relax_T = 1.0;
    // solve_gssor(nx, ny, aP, aE, aW, aN, aS, b, T, Tpad, Tpnew, max_iter, tol, relax_T);
    // long long timeTaken = solve_gssor_wave(nx, ny, aP, aE, aW, aN, aS, b, T, Tpad, Tpnew, max_iter, tol, relax_T, num_threads);
    long long timeTaken = solve_gssor_wave();

    // ---
    // printf(" > Done solving for T ------------- \n");

    printf("\nNumber of iterations: %d\n", num_iters_taken);
    printf("Final error: %9.5e\n", rel_err);
    printf("Time taken: %lld ms\n\n", timeTaken);

    get_exact_soln(nx, ny, x, y, Tex);
    output_soln(nx, ny, 0, x, y, T, Tex);

    // free memory
    // ----1D arrays ---
    free(y);
    free(dy);
    free(x);
    free(dx);
    // --- Done 1D arrays ---

    // ----2D arrays ---
    for (i = 0; i < nx; i++)
        free(T[i]);
    free(T);
    for (i = 0; i < nx; i++)
        free(Tex[i]);
    free(Tex);
    for (i = 0; i < nx; i++)
        free(b[i]);
    free(b);
    for (i = 0; i < nx; i++)
        free(aP[i]);
    free(aP);
    for (i = 0; i < nx; i++)
        free(aE[i]);
    free(aE);
    for (i = 0; i < nx; i++)
        free(aW[i]);
    free(aW);
    for (i = 0; i < nx; i++)
        free(aN[i]);
    free(aN);
    for (i = 0; i < nx; i++)
        free(aS[i]);
    free(aS);
    for (i = 0; i < nx; i++)
        free(kdiff[i]);
    free(kdiff);
    for (i = 0; i < ny; i++)
        free(bcleft[i]);
    free(bcleft);
    for (i = 0; i < ny; i++)
        free(bcrght[i]);
    free(bcrght);
    for (i = 0; i < nx; i++)
        free(bctop[i]);
    free(bctop);
    for (i = 0; i < nx; i++)
        free(bcbot[i]);
    free(bcbot);
    for (i = 0; i < nx + 2; i++)
        free(Tpnew[i]);
    free(Tpnew);
    free(thread_ids);
    free(thread_errs);
    free(thread_tpad_max);
    free(thread_tpnew_max);

    pthread_barrier_destroy(&bar_calc);
    pthread_barrier_destroy(&bar_upd);

    // --- Done 2D arrays ---
    // printf("\n > Done freeing up memory --------- \n");
    return 0;
}