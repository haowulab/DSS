/* utility functions for filtering */
#include "DSS.h"


/* engine function for moving avearge.
   Note this function uses cumsum results, so it's very fast. 
   flag=0 means to do moving sum, flag=1 means moving average
   Compared to "filter" function in R this can *only* smooth with rectangular kernel. */
void windowFilter_engine(double *x, int* nitem, int npos, int flag, double* result_ptr) {
  int ipos, j;
  double *sx;
  SEXP result;

  /* first compute cumsums. This will greatly speed up the calculation */
  sx = (double *)R_alloc(npos+1, sizeof(double));
  sx[0]=0;
  for(j=0; j<npos; j++)
    sx[j+1] = sx[j] + x[j];

  /* Now perform moving average. Be careful at boundaries. */
  for(ipos=0; ipos<npos; ipos++) {
    if(flag==0)
      result_ptr[ipos]=(sx[ipos+1+nitem[ipos+npos]]-sx[ipos-nitem[ipos]]);
    if(flag==1)
      result_ptr[ipos]=(sx[ipos+1+nitem[ipos+npos]]-sx[ipos-nitem[ipos]]) / (nitem[ipos+npos]+nitem[ipos]+1);
  }
}

/* wrapper function to do moving average within a fixed window size */
SEXP windowFilter(SEXP x, SEXP pos, SEXP ws, SEXP R_flag) {
  SEXP result;
  int npos = length(pos);
  int *nitem=(int *)R_alloc(npos*2, sizeof(int));
  int flag=*INTEGER(R_flag);
  nitem_bin(pos, *INTEGER(ws), nitem);

  PROTECT(result = allocVector(REALSXP, npos));
  double *result_ptr = REAL(result);
  windowFilter_engine(REAL(x), nitem, npos, flag, result_ptr);

  UNPROTECT(1);
  return(result);
}

