#include "mygemm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

#include "util.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

void dgemm0(const double* A, const double* B, double* C, const int n)
{
	int i=0;
	int j=0;
	int k=0;
	
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			for (k=0; k<n; k++)
				C[i*n+j] += A[i*n+k] * B[k*n+j];
//	printf("mygemm c = %lf\n", C[0]);

//	int curr_dim = n;
//            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, curr_dim, curr_dim, curr_dim, 1., A, curr_dim, B, curr_dim, 1., C, curr_dim);

}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
	int i = 0;
	int j = 0;
	int k = 0;

	for (i = 0; i < n; i++)
		for (j = 0; j <n ; j++) 
		{
			register double r = C[i*n+j] ;
			for (k = 0; k < n; k++)
				r += A[i*n+k] * B[k*n+j];
			C[i*n+j] = r;
		}		
}

void dgemm2(const double *A, const double *B, double *C, const int n) 
{
	int i, j, k;
	for (i = 0; i < n; i+=2)
		for (j = 0; j < n; j+=2)
			for (k = 0; k<n; k+=2)
			{
			C[i*n + j] 		= A[i*n + k] * B[k*n + j] 	+ A[i*n + k+1]*B[(k+1)*n + j] 		+ C[i*n + j];
			C[(i+1)*n + j] 		= A[(i+1)*n + k]*B[k*n + j] 	+ A[(i+1)*n + k+1]*B[(k+1) * n + j] 	+ C[(i+1)*n + j];
			C[i*n + (j+1)]  	= A[i*n + k]*B[k*n + (j+1)] 	+ A[i*n + k+1]*B[(k+1) * n + (j+1)] 	+ C[i*n + (j+1)];
			C[(i+1)*n + (j+1)] 	= A[(i+1)*n + k]*B[k*n + (j+1)] + A[(i+1)*n + k+1] * B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
			}

}

void dgemm3(const double *A, const double *B, double *C, const int n) 
{

	int i , j , k;
	for(i = 0; i < n; i += 2)
       		for(j = 0; j < n; j += 2)  
		{ 
            	register double C00 = C[i*n+j]; register double C01 = C[i*n+j+1];  register double C10 = C[(i+1)*n+j]; register double C11 = C[(i+1)*n+j+1];

            	for(k = 0; k < n; k += 2) 
		{
                	register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
	                register double a00 = A[ta]; register double a01 = A[ta+1]; register double a10 = A[tta]; register double a11 = A[tta+1];
        	        register double b00 = B[tb]; register double b01 = B[tb+1]; register double b10 = B[ttb]; register double b11 = B[ttb+1];
	                C00 += a00*b00 + a01*b10;
        	        C01 += a00*b01 + a01*b11;
                	C10 += a10*b00 + a11*b10;
	                C11 += a10*b01 + a11*b11;
        	}

             	C[i*n+j] = C00;
	        C[i*n+j+1] = C01;
        	C[(i+1)*n+j] = C10;
             	C[(i+1)*n+j+1] = C11;
       }

}

void ijk(const double *A, const double *B, double *C, const int n) 
{

}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jik(const double *A, const double *B, double *C, const int n) 
{

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kij(const double *A, const double *B, double *C, const int n) 
{

}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{

}


void ikj(const double *A, const double *B, double *C, const int n) 
{

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void jki(const double *A, const double *B, double *C, const int n) 
{

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kji(const double *A, const double *B, double *C, const int n) 
{

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}
