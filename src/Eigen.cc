/* eigen.c */



#include "TransMat.h"
#include "Eigen.h"

#define HAVE_LAPACK

extern "C" {

SEXP lambda(SEXP x)
{
    int i, n, lwork, info, *xdims;
    double *work, *wR, *wI, *left, *right, *xvals, tmp, maxval;
    char jobVL[1], jobVR[1];
    SEXP ret;

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1])
	error("x must be a square numeric matrix");

    xvals = (double *) R_alloc(n * n, sizeof(double));
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), (size_t) (n * n));

    for (i=0; i<(n*n); i++)
      Rprintf("%i %f \n",i,REAL(x)[i]);

    //    vectors = !ov;
    jobVL[0] = jobVR[0] = 'N';
    left = right = (double *) 0;
    wR = (double *) R_alloc(n, sizeof(double));
    wI = (double *) R_alloc(n, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
#ifdef HAVE_LAPACK
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, &tmp, &lwork, &info);
#else
    F77_CALL(rgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, &tmp, &lwork, &info);
#endif
    if (info != 0)
	error("error code %d from Lapack routine dgeev", info);
    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));
#ifdef HAVE_LAPACK
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, work, &lwork, &info);
#else
    F77_CALL(rgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, work, &lwork, &info);
#endif
    if (info != 0)
	error("error code %d from Lapack routine dgeev", info);

    ret = PROTECT(allocVector(REALSXP, 1));
    maxval=-1000000.0;
    for (i = 0; i < n; i++)
      {
	if (wI[i]==0.0)
	  {
	    if (wR[i]>maxval) maxval=wR[i];
	  }

      }

    REAL(ret)[0] = maxval;

    UNPROTECT(1);
    return ret;
}


} //extern "C"



