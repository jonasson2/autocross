#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declare the Fortran subroutines */
void F77_NAME(p3_subroutine)(int *n, double *time, double *x, double *y,
                             double *alpha, double *corr, double *ci,
                             double *taux, double *tauy);

void F77_NAME(rx_subroutine)(int *nx, int *ny, int *nout,
                             double *tx, double *x, double *ty, double *y,
                             int *cfg_nsim, double *cfg_ofac, double *cfg_hifac,
                             int *cfg_n50, double *cfg_alpha, int *cfg_iwin,
                             double *rhox, double *rhoy, double *taux, double *tauy,
                             double *dof, double *db6, double *false_alarm,
                             double *faccritx, double *faccrity, 
                             double *alphacritx, double *alphacrity,
                             double *data_x, double *data_y,
                             double *data_xy, double *data_cxy, double *data_phxy);

void F77_NAME(rx_setdim)(int *nx, int *ny,
                         double *tx, double *ty,
                         double *cfg_ofac, double *cfg_hifac,
                         int *cfg_n50, int *nout);

/* Register the Fortran routines */
static const R_FortranMethodDef FortranEntries[] = {
    {"p3_subroutine", (DL_FUNC) &F77_NAME(p3_subroutine), 9},
    {"rx_subroutine", (DL_FUNC) &F77_NAME(rx_subroutine), 29},
    {"rx_setdim", (DL_FUNC) &F77_NAME(rx_setdim), 8},
    {NULL, NULL, 0}
};

/* Initialize the package */
void R_init_autocross(DllInfo *dll) {
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
