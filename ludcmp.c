
#include "decs.h"
#define TINY 1.0e-30;

int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d)
{
  int i, imax, j, k;
  FTYPE big, dum, sum, temp;
  static FTYPE *vv;
  static int firstc = 1;

  imax = 0;
  if (firstc) {
    firstc = 0;
    vv = dvector(1, n);
  }
  *d = 1.0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((temp = fabs(a[i][j])) > big)
	big = temp;
    if (big == 0.0) {
      dualfprintf(fail_file, "sing. matr. in ludcmp: %d %d %d\n", icurr,jcurr,kcurr);
      for (i = 1; i <= n; i++) for(j=1;j<=n;j++){
	dualfprintf(fail_file,"a[%d][%d]=%21.15g\n",i,j,a[i][j]);
	if(!finite(a[i][j])){
	  dualfprintf(fail_file,"i=%d j=%d is not finite in ludcmp: icurr=%d jcurr=%d\n",i,j,icurr,jcurr);
	}
      }
      return(1);
      /* 
         nrerror("Singular matrix in routine ludcmp"); */
    }
    vv[i] = 1.0 / big;
  }
  for (j = 1; j <= n; j++) {
    for (i = 1; i < j; i++) {
      sum = a[i][j];
      for (k = 1; k < i; k++)
	sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++) {
      sum = a[i][j];
      for (k = 1; k < j; k++)
	sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k = 1; k <= n; k++) {
	dum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = TINY;
    if (j != n) {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i <= n; i++)
	a[i][j] *= dum;
    }
  }
  return(0);
}

#undef TINY
