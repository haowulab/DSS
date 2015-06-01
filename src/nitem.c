/* find number of items in each bin.
   result is a two column matrix for number of items
   before and after the current position. */

#include "DSS.h"

void nitem_bin(SEXP pos_sexp, int ws, int* result_ptr) {
  int ipos, j, n;
  SEXP result;
  int npos = length(pos_sexp);
  int *pos = INTEGER(pos_sexp);

  for(ipos=0; ipos<npos; ipos++) {
    /* left side */
    n = 0;
    for(j=ipos-1; j>=0; j--) {
      if(pos[ipos]-pos[j]>ws/2) {
	n = ipos - j - 1; 
	break;
      }
    }
    if(j == -1) /* begining */
      n = ipos; 
    result_ptr[ipos] = n;

    /* right side */
    n = 0;
    for(j=ipos+1; j<npos; j++) {
      if(pos[j]-pos[ipos]>ws/2) {
        n = j-ipos-1;
	break;
      } 
    }
    if(j == npos) /* end */
      n = npos - ipos - 1;
    result_ptr[ipos+npos] = n;
  }
}
