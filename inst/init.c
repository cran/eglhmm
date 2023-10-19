#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(esttpm)(int *ns, int *n, int *k, double *tpm,
                             int *mixture, double *wrk);

extern void F77_NAME(getgl)(int *ndistr, double *fy, double *zeta,
                            int *nz, double *phimat, double *y,
                            int *ymiss, double *dmat, int *nind,
                            int *n, double *tpm, double *xispd,
                            double *sigma, int *size, int *nbot, int *ntop,
                            double *d1pi,
                            int *kstate, int *npar, int *nphi,
                            int *nyv, int *nxc,
                            double *d1p, double *d1f,
                            double *d1zeta, double *d1u,
                            double *d1a, double *d1b,
                            double *alpha, double *alphw, double *a,
                            double *aw, double *xlc);

extern void F77_NAME(gethgl)(int *ndistr, double *fy, double *zeta,
                             int *nz,
                             double *phimat, double *y,
                             int *ymiss, double *dmat, int *nind, int *n,
                             double *tpm, double *xispd, double *sigma,
                             int *size, int *nbot, int *ntop, double *d1pi,
                             double *d2pi, int *kstate,
                             int *npar, int *nphi, int *nyv, int *nxc,
                             double *d1p, double *d2p,
                             double *d1f, double *d2f, double *d1zeta,
                             double *d2zeta, double *d1u, double *d2u,
                             double *d1a, double *d1b, double * d2aa,
                             double *d2ab, double *d2bb,
                             double *alpha, double *alphw, double *a, double *b,
                             double *aw, double *bw, double *xlc, double *hess);

extern void F77_NAME(getll)(double *fy, double *tpm,
                           double *xispd, int *kstate, int *n, double *alpha,
                           double *alphw, double *xlc);

extern void F77_NAME(recurse)(double *fy, double *xispd, double *tpm,
                              double *epsilon, int *kstate, int *nc,
                              double *wrk, double *xlc, double *alpha,
                              double *beta, double *gamma, double *xi,
                              double *xisum);

static const R_FortranMethodDef FortranEntries[] = {
    {"esttpm",  (DL_FUNC) &F77_NAME(esttpm),   6},
    {"getgl",   (DL_FUNC) &F77_NAME(getgl),   33},
    {"gethgl",  (DL_FUNC) &F77_NAME(gethgl),  44},
    {"getll",   (DL_FUNC) &F77_NAME(getll),    8},
    {"recurse", (DL_FUNC) &F77_NAME(recurse), 13},
    {NULL, NULL, 0}
};

void R_init_eglhmm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
