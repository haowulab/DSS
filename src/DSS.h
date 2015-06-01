#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#define BINSIZE 50

SEXP compute_var_smooth(SEXP vars, SEXP delta, SEXP n1, SEXP pos, SEXP ws, SEXP rho);
SEXP windowFilter(SEXP x, SEXP pos, SEXP ws, SEXP R_flag);
void compute_var_smooth_engine(double *vars, double *delta, double *n1, int *nitem,
			   double *result_ptr, int npos, double rho, int *pos);
void windowFilter_engine(double *x, int* nitem, int npos, int flag, double* result_ptr);
void nitem_bin(SEXP pos_sexp, int ws, int* result_ptr);

