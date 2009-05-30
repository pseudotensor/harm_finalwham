#include "decs.h"

void nrerror(char error_text[])
{
  dualfprintf(fail_file, "Numerical Recipes run-time error...\n");
  dualfprintf(fail_file, "%s\n", error_text);
  dualfprintf(fail_file, "...now exiting to system...\n");
  exit(1);
}



float *vector(long nl, long nh)
{
  float *v;

  v = (float *) malloc((unsigned) (nh - nl + 1) * sizeof(float));
  if (!v)
    nrerror("allocation failure in vector()");
  return v - nl;
}

int *ivector(long nl, long nh)
{
  int *v;

  v = (int *) malloc((unsigned) (nh - nl + 1) * sizeof(int));
  if (!v)
    nrerror("allocation failure in ivector()");
  return v - nl;
}

FTYPE *dvector(long nl, long nh)
{
  FTYPE *v;

  v = (FTYPE *) malloc((unsigned) (nh - nl + 1) * sizeof(FTYPE));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v - nl;
}



float **matrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  float **m;

  m = (float **) malloc((unsigned) (nrh - nrl + 1) * sizeof(float *));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (float *) malloc((unsigned) (nch - ncl + 1) * sizeof(float));
    if (!m[i])
      nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  FTYPE **m;

  m = (FTYPE **) malloc((unsigned) (nrh - nrl + 1) * sizeof(FTYPE *));
  if (!m)
    nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] =
	(FTYPE *) malloc((unsigned) (nch - ncl + 1) * sizeof(FTYPE));
    if (!m[i])
      nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  int **m;

  m = (int **) malloc((unsigned) (nrh - nrl + 1) * sizeof(int *));
  if (!m)
    nrerror("allocation failure 1 in imatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (int *) malloc((unsigned) (nch - ncl + 1) * sizeof(int));
    if (!m[i])
      nrerror("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}



float **submatrix(float **a, long oldrl, long oldrh, long oldcl,
		  long oldch, long newrl, long newcl)
{
  long i, j;
  float **m;

  m = (float **) malloc((unsigned) (oldrh - oldrl + 1) *
			sizeof(float *));
  if (!m)
    nrerror("allocation failure in submatrix()");
  m -= newrl;

  for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
    m[j] = a[i] + oldcl - newcl;

  return m;
}



void free_vector(float *v, long nl, long nh)
{
  free((char *) (v + nl));
}

void free_ivector(int *v, long nl, long nh)
{
  free((char *) (v + nl));
}

void free_dvector(FTYPE *v, long nl, long nh)
{
  free((char *) (v + nl));
}



void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}

void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  long i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}



void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
{
  free((char *) (b + nrl));
}



float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
{
  long i, j, nrow, ncol;
  float **m;

  nrow = nrh - nrl + 1;
  ncol = nch - ncl + 1;
  m = (float **) malloc((unsigned) (nrow) * sizeof(float *));
  if (!m)
    nrerror("allocation failure in convert_matrix()");
  m -= nrl;
  for (i = 0, j = nrl; i <= nrow - 1; i++, j++)
    m[j] = a + ncol * i - ncl;
  return m;
}



void free_convert_matrix(float **b, long nrl, long nrh, long ncl,
			 long nch)
{
  free((char *) (b + nrl));
}
