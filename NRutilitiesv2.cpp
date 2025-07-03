/*
  NRutilities.cpp

  Selected utility functions from Numerical Recipes.
  These are in the public domain
*/

#include "NRutilities.h"
#include <cstdio>
#include <cstddef>
#include <cstdlib>

#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	std::fprintf(stderr, "Numerical Recipes run-time error...\n");
	std::fprintf(stderr, "%s\n", error_text);
	std::fprintf(stderr, "...now exiting to system...\n");
	std::exit(1);
}

/*******************************************************************/

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v = static_cast<unsigned char *>(std::malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned char))));
	if (!v) nrerror("allocation failure in cvector()");
	return v - nl + NR_END;
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v = static_cast<float *>(std::malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float))));
	if (!v) nrerror("allocation failure in vector()");
	return v - nl + NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v = static_cast<int *>(std::malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int))));
	if (!v) nrerror("allocation failure in ivector()");
	return v - nl + NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v = static_cast<double *>(std::malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double))));
	if (!v) nrerror("allocation failure in dvector()");
	return v - nl + NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	/* allocate pointers to rows */
	m = static_cast<float **>(std::malloc((size_t)((nrow + NR_END) * sizeof(float*))));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = static_cast<float *>(std::malloc((size_t)((nrow*ncol + NR_END) * sizeof(float))));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;

	/* allocate pointers to rows */
	m = static_cast<double **>(std::malloc((size_t)((nrow + NR_END) * sizeof(double*))));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = static_cast<double *>(std::malloc((size_t)((nrow*ncol + NR_END) * sizeof(double))));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int **m;

	/* allocate pointers to rows */
	m = static_cast<int **>(std::malloc((size_t)((nrow + NR_END) * sizeof(int*))));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl] = static_cast<int *>(std::malloc((size_t)((nrow*ncol + NR_END) * sizeof(int))));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] 
   This function is not in NRutils, but it uses the same structure */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	char **m;

	/* allocate pointers to rows */
	m = static_cast<char **>(std::malloc((size_t)((nrow + NR_END) * sizeof(char*))));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl] = static_cast<char *>(std::malloc((size_t)((nrow*ncol + NR_END) * sizeof(char))));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

/*******************************************************************/

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	std::free(v + nl - NR_END);
}


void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	std::free(v + nl - NR_END);
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	std::free(v + nl - NR_END);
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	std::free(v + nl - NR_END);
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	std::free(m[nrl] + ncl - NR_END);
	std::free(m + nrl - NR_END);
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	std::free(m[nrl] + ncl - NR_END);
	std::free(m + nrl - NR_END);
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	std::free(m[nrl] + ncl - NR_END);
	std::free(m + nrl - NR_END);
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a char matrix allocated by cmatrix() */
/* BUG FIX: The original file incorrectly listed the first argument as int **m */
{
	std::free(m[nrl] + ncl - NR_END);
	std::free(m + nrl - NR_END);
}

/*******************************************************************/