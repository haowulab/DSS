/* function to compute variances when there's smoothing. 
   This works for one chr only */
#include "DSS.h"

/* wrapper function */
SEXP compute_var_smooth(SEXP vars, SEXP delta, SEXP n1, SEXP pos, SEXP ws, SEXP rho) {
	SEXP result;
	int npos = length(pos);

	/* compute number of CG sites in each window */
	int *nitem=(int *)R_alloc(npos*2, sizeof(int));
	nitem_bin(pos, *INTEGER(ws), nitem);

	PROTECT(result = allocVector(REALSXP, npos));
	double *result_ptr = REAL(result);
	
	/* compute variances */
	compute_var_smooth_engine(REAL(vars), REAL(delta), REAL(n1), nitem, result_ptr, 
							  npos, *REAL(rho), INTEGER(pos));
	UNPROTECT(1);
	return(result);

}

/* engine function for computing variances */
void compute_var_smooth_engine(double *vars, double *delta, double *n1, int *nitem, 
							   double *result_ptr, int npos, double rho, int *pos) {
	int ipos, istart, iend, j, k, nCG_bin, nlag;
	double thisvar, dtmp; 

	/* loop on CG site */
	for(ipos=0; ipos<npos; ipos++) {
		istart = ipos-nitem[ipos];
		iend = ipos+nitem[ipos+npos];
		nCG_bin = iend - istart + 1;
		thisvar = 0;

		if(nCG_bin <= 1) /* only one CG site in this window */
			thisvar = vars[ipos];
		else { /* multiple CG sites in the window, need to compute var/cov matrix */
			for( j=istart; j<iend; j++ ) {
				for( k=j; k<iend; k++) {
					if(j == k) 
						thisvar = thisvar + vars[ipos];
					else { /* different position */
						/* distances. This is used to compute correlation. */
						nlag = abs(pos[k] - pos[j]) / BINSIZE;
						dtmp = n1[j]*n1[k]*pow(rho, nlag) * (delta[j] * delta[k]);
						thisvar = thisvar + dtmp*2;
					}
				}
			}
		}
		result_ptr[ipos ] = thisvar;
	}
}
