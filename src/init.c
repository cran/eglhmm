#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(esttpm)(int *ns, int *n, int *k, double *tpm,
                             int *mixture, double *wrk);

extern void F77_NAME(gethgl)(int *nd, int *ndistr, double *fy, double *gmu,
                             double *sd, double *lambda, double *p,
                             double *ashp, double *bshp, double *phimat,
                             double *y, int *ymiss, double *tdm, int *nind, int *n,
                             double *tpm, double *xispd, int *size, int *nbot,
                             int *ntop, double *d1pi, double *d2pi, int *kstate,
                             int *npar, int *npt, int *nyv, int *nxc, double *d1p,
                             double *d2p, double *d1f, double *d2f, double *d1a,
                             double *d1b, double *d2aa, double *d2ab,
                             double *d2bb, double *alpha, double *alphw,
                             double *a, double *b, double *aw, double *bw,
                             double *xlc, double *hess);

extern void F77_NAME(recurse)(double *fy, double *xispd, double *tpm,
                              double *epsilon, int *kstate, int *nc,
                              double *wrk, double *xlc, double *alpha,
                              double *beta, double *gamma, double *xi,
                              double *xisum, int *level);

static const R_FortranMethodDef FortranEntries[] = {
    {"esttpm",  (DL_FUNC) &F77_NAME(esttpm),   6},
    {"gethgl",  (DL_FUNC) &F77_NAME(gethgl),  44},
    {"recurse", (DL_FUNC) &F77_NAME(recurse), 14},
    {NULL, NULL, 0}
};

void R_init_eglhmm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
