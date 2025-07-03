/*
matrix.cpp   code for some matrix procedures

This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by 
the Free Software Foundation.

This file is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public 
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>   

Copyright 2009 
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa

*/

#include <cmath>
#include "matrix.h" 
#include "NRutilities.h"
#include "ERutilities.h"

/*--------------------------------------------------------------------------
	MatrixMultiVector0
	
	functionality:

	Computes a nrow x ncol matrix multipled by a ncol vector, the matrix and
	vector are 0 offset.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow        number of rows
		ncol        number of columns
		m           nrow x ncol matrix
		v           ncol vector

 
	output:
	    r           resulting vector of nrow elements
--------------------------------------------------------------------------*/
void MatrixMultiVector0(int nrow, int ncol, double **m, double *v, double *r)
{
	int i, j;

	for (i=0; i<nrow; i++) {
		r[i] = 0;
		for (j=0; j<ncol; j++) 
			r[i] += m[i][j] * v[j];
	}
}


/*--------------------------------------------------------------------------
	MatrixMultiVector1
	
	functionality:

	Computes a nrow x ncol matrix multipled by a ncol vector, the matrix and
	vector are 1 offset.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow        number of rows
		ncol        number of columns
		m           nrow x ncol matrix
		v           ncol vector

 
	output:
	    r           resulting vector of nrow elements
--------------------------------------------------------------------------*/
void MatrixMultiVector1(int nrow, int ncol, double **m, double *v, double *r)
{
	int i, j;

	for (i=1; i<=nrow; i++) {
		r[i] = 0;
		for (j=1; j<=ncol; j++) 
			r[i] += m[i][j] * v[j];
	}
}


/*--------------------------------------------------------------------------
	VectorMultiMatrix0
	
	functionality:

	Computes a nrow vector  multipled by a nrow x ncol matrix, the matrix and
	vector are 0 offset
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow        number of rows
		ncol        number of columns
		m           nrow x ncol matrix
		v           nrow vector

 
	output:
	    r           resulting vector of ncol elements
--------------------------------------------------------------------------*/
void VectorMultiMatrix0(int nrow, int ncol, double *v, double **m, double *r)
{
	int i, j;

	for (i=0; i<ncol; i++) {
		r[i] = 0;
		for (j=0; j<nrow; j++) r[i] += v[j] * m[j][i];
	}
}

/*--------------------------------------------------------------------------
	VectorMultiMatrix1
	
	functionality:

	Computes a nrow vector  multipled by a nrow x ncol matrix, the matrix and
	vector are 1 offset
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow        number of rows
		ncol        number of columns
		m           nrow x ncol matrix
		v           nrow vector

 
	output:
	    r           resulting vector of ncol elements
--------------------------------------------------------------------------*/
void VectorMultiMatrix1(int nrow, int ncol, double *v, double **m, double *r)
{
	int i, j;

	for (i=1; i<=ncol; i++) {
		r[i] = 0;
		for (j=1; j<=nrow; j++) r[i] += v[j] * m[j][i];
	}
}


/*--------------------------------------------------------------------------
	MatrixMultiMatrix0
	
	functionality:

	Computes m (nrow1 x ncol1row2 matrix) multipled by n ( ncol1row2 x ncol2 matrix), 
	the matrix and vector are 0 offset.
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow1       number of rows for m
		ncol1row2   number of columns for m and number of rows for n
		ncol2       number of columnss for n
		m           nrow1 x ncol matrix
		n           ncol x ncol2 matrix

 
	output:
	    r           resulting matrix of nrow1 x ncol2
--------------------------------------------------------------------------*/
void MatrixMultiMatrix0(int nrow1, int ncol1row2, int ncol2, double **m, double **n, double **r)
{
	int i, j, k;

	for (i=0; i<nrow1; i++) 
		for (j=0; j<ncol2; j++) r[i][j] = 0;
	
	for (i=0; i<nrow1; i++) 
		for (j=0; j<ncol2; j++) 
			for (k=0; k<ncol1row2; k++)  
				r[i][j] += m[i][k] * n[k][j];
	
}


/*--------------------------------------------------------------------------
	MatrixMultiMatrix1
	
	functionality:

	Computes m (nrow1 x ncol1row2 matrix) multipled by n ( ncol1row2 x ncol2 matrix),
	, the matrix and vector are 1 offset.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow1       number of rows for m
		ncol1row2   number of columns for m and number of rows for n
		ncol2       number of columnss for n
		m           nrow1 x ncol matrix
		n           ncol x ncol2 matrix

 
	output:
	    r           resulting matrix of nrow1 x ncol2
--------------------------------------------------------------------------*/
void MatrixMultiMatrix1(int nrow1, int ncol1row2, int ncol2, double **m, double **n, double **r)
{
	int i, j, k;

	for (i=1; i<=nrow1; i++) 
		for (j=1; j<=ncol2; j++) r[i][j] = 0;
	
	for (i=1; i<=nrow1; i++) 
		for (j=1; j<=ncol2; j++) 
			for (k=1; k<=ncol1row2; k++)  
				r[i][j] += m[i][k] * n[k][j];
	
}

/*--------------------------------------------------------------------------
	VectorNormSq
	
	functionality:

	Computes the square of the norm of a vector.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		ncat       Number of discrete score categories 
		v          vector with ncat elements
 
 
	output:
 		return the norm of the vector
--------------------------------------------------------------------------*/
double VectorNorm0(int ncat, double *v)
 {
	 int i;
	 double sum=0;

	 for (i=0; i<ncat; i++)	
		 sum += v[i] * v[i];

	 return std::sqrt(sum);
 }

 /*******************************************************************************
PTransposeMultiP1

Functionality: Compute the transpose of a matrix multiplied by itself

Author: Tianyou Wang 

Date:   1/11/06

Input:
       nr       number of rows of P 
	   nc       number of columns of P 
	   P        the input matrix, index from 1 to nr, 1 to nc.

Output:
       PtP      the product matrix

********************************************************************************/

void PTransposeMultiP1(int nr, int nc, double **P, double **PtP)
{
	int i, j, k;
	double sum;

	for (i=1; i<=nc; i++) {
		for (j=1; j<=nc; j++) {
			sum = 0;
			for (k=1; k<=nr; k++) {
				sum += P[k][i] * P[k][j];
			}
			PtP[i][j] = sum;
		}
	}
}

 /*******************************************************************************
PTransposeMultiP0

Functionality: Compute the transpose of a matrix multiplied by itself

Author: Tianyou Wang 

Date:   1/11/06

Input:
       nr       number of rows of P 
	   nc       number of columns of P 
	   P        the input matrix, index from 1 to nr, 1 to nc.

Output:
       PtP      the product matrix

********************************************************************************/

void PTransposeMultiP0(int nr, int nc, double **P, double **PtP)
{
	int i, j, k;
	double sum;

	for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
			sum = 0;
			for (k=0; k<nr; k++) {
				sum += P[k][i] * P[k][j];
			}
			PtP[i][j] = sum;
		}
	}
}
				

/*********************************************************************
Determinant

Functionality: compute the determinant of a square matrix

Author: Tianyou Wang

Date:   4/19/05

Input:
       a     matrix for which the determinant is computed
	   n     the dimension of the matrix

Output:
       the function return the determinant of matrix a

*********************************************************************/
double Determinant(double **a, int n)
{
	int i, j, *indx;
	double **aa;
	double d = 1;

	/* because ludcmp takes a 1-ofset matrix but the input matrix is 
	   0-ofset, the following code are to copy a into aa */
	aa = dmatrix(1, n, 1, n);
	indx = ivector(1, n);

	for(i=1; i<=n; i++)
		for(j=1; j<=n; j++)
			aa[i][j] = a[i-1][j-1];
	
	er_ludcmp(aa, n, indx, &d); 
	d=1;
	
	for(j=1; j<=n; j++) d *= aa[j][j];

    // BUG FIX: Free allocated memory
    free_dmatrix(aa, 1, n, 1, n);
    free_ivector(indx, 1, n);

	return d;
}

/*********************************************************************
MatrixInverse1

Functionality: compute the inverse of a square matrix

Author: Tianyou Wang

Date:   1/12/06

Input:
       m     matrix for which the determinant is computed
	   n     the dimension of the matrix

Output:
       m     the inverse matrix of m
*********************************************************************/
/* ===== Before version 1.0 update ===== 
void MatrixInverse1(int n, double **m)
{
	int i, j;
	double **Identity;

	Identity = dmatrix(1, n, 1, n);
	for(i=1; i<=n; i++) {
		for(j=1; j<=n; j++) {
			Identity[i][j] = 0;
		}
		Identity[i][i] = 1;
	}

	gaussj(m, n, Identity , n);  

}
*/

/*********************************************************************
MatrixInverse0

Functionality: compute the inverse of a square matrix

Author: Tianyou Wang

Date:   1/12/06

Input:
       m     matrix for which the determinant is computed
	   n     the dimension of the matrix

Output:
       m     the inverse matrix of m
*********************************************************************/
/* ===== Before version 1.0 update ===== 
void MatrixInverse0(int n, double **m)
{
	int i, j;
	double **m1;

	m1 = dmatrix(1, n, 1, n);
	for(i=1; i<=n; i++) {
		for(j=1; j<=n; j++) {
			m1[i][j] = m[i-1][j-1];
		}
	}

	MatrixInverse1(n, m1);  

	for(i=1; i<=n; i++) {
		for(j=1; j<=n; j++) {
			m[i-1][j-1] = m1[i][j];
		}
	}
}
*/


 /*******************************************************************************
PtAP0

Functionality: Compute the transpose of a P multiplied by square matrix A 
and then multiplied by P

Author: Tianyou Wang 

Date:   1/11/06

Input:
       nr       number of rows of P 
	   nc       number of columns of P (usually nc < nr)
	   P        the input matrix, index from 1 to nr, 1 to nc.
	   A        matrix (nr x nr) 

Output:
       PtAP      the product matrix (nc x nc)

********************************************************************************/

void PtAP0(int nr, int nc, double **P, double **A, double **PtAP)
{
	int i, j, k, l;
	double sum;
	double *temp;

	temp = dvector(0, nr-1); // Corrected size

	for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
			sum = 0;
			for (l=0; l<nr; l++) {
				temp[l] = 0;
				for (k=0; k<nr; k++) {
					temp[l] += P[k][i] * A[k][l];
				}
				sum += temp[l] * P[l][j];
			}
			PtAP[i][j] = sum;
		}
	}
    
    // BUG FIX: Free allocated memory
    free_dvector(temp, 0, nr-1);
}
				
 /*******************************************************************************
vtAv0

Functionality: Compute the transpose of a vector v multiplied by square matrix A 
and then multiplied by v

Author: Tianyou Wang 

Date:   1/11/06

Input:
       n		number of elements in v 
	   nc       number of columns of P (usually nc < nr)
	   v        the input vector.
	   A        matrix (n x n) 

Output:
       vtAv      the product matrix (nc x nc)

********************************************************************************/

double vtAv0(int n, double *v, double **A)
{
	int k, l;
	double sum;
	double *temp;

	temp = dvector(0, n-1); // Corrected size

	sum = 0;
	for (l=0; l<n; l++) {
		temp[l] = 0;
		for (k=0; k<n; k++) {
			temp[l] += v[k] * A[k][l];
		}
		sum += temp[l] * v[l];
	}

    // BUG FIX: Free allocated memory
    free_dvector(temp, 0, n-1);

	return sum;
}