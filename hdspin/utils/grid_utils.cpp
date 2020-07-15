/* Auxiliary files grid utilities.
 *
 * Matthew Carbone, Columbia University 2020
 */

#include <map>
#include <vector>
#include <math.h>

#include "grid_utils.h"

void fill_pyLogspace(double *arr, double a, double b, int n_points,
    const double base)
{
    for (int ii = 0; ii < n_points; ii++)
    {
        arr[ii] = pow(base, ii * (b - a) / (n_points - 1));
    }
}

/* Get's the first pi_basin/config grid. */
void fill_pi_grid_1(double *arr, const double nMC, const double dw,
    const int n_samp)
{
    const double tw_max = nMC / (dw + 1.0);
    fill_pyLogspace(arr, 0, log10(tw_max), n_samp, 10.0);
    arr[0] = 0.0;  // Start with tw(0) = 0
}

/* Get's the second pi_basin/config grid from the first one. */
void  fill_pi_grid_2(double *arr, const double *pi_grid_1,
    const double dw, const int n_samp)
{
    for (int ii=0; ii<n_samp; ii++){arr[ii] = pi_grid_1[ii] * (dw + 1.0);}
}






