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
int BLOCK = 256;
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
	/*int i, j, k;
	for (i = 0; i < n; i+=2)
		for (j = 0; j < n; j+=2)
			for (k = 0; k<n; k+=2)
			{
			C[i*n + j] 		= A[i*n + k] * B[k*n + j] 	+ A[i*n + k+1]*B[(k+1)*n + j] 		+ C[i*n + j];
			C[(i+1)*n + j] 		= A[(i+1)*n + k]*B[k*n + j] 	+ A[(i+1)*n + k+1]*B[(k+1) * n + j] 	+ C[(i+1)*n + j];
			C[i*n + (j+1)]  	= A[i*n + k]*B[k*n + (j+1)] 	+ A[i*n + k+1]*B[(k+1) * n + (j+1)] 	+ C[i*n + (j+1)];
			C[(i+1)*n + (j+1)] 	= A[(i+1)*n + k]*B[k*n + (j+1)] + A[(i+1)*n + k+1] * B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
			}*/
        int i, j , k;
        for(i = 0; i < n; i += 2)
                for(j = 0; j < n; j += 2)  {
                        register int t = i*n+j; register int tt = t+n;
                        register double c00 = C[t]; register double c01 = C[t+1];  register double c10 = C[tt]; register double c11 = C[tt+1];
                        for(k = 0; k < n; k += 2) {
                                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                                register double a00 = A[ta]; register double a01 = A[ta+1]; register double a10 = A[tta]; register double a11 = A[tta+1];
                                register double b00 = B[tb]; register double b01 = B[tb+1]; register double b10 = B[ttb]; register double b11 = B[ttb+1];
                                c00 += a00*b00 + a01*b10;
                                c01 += a00*b01 + a01*b11;
                                c10 += a10*b00 + a11*b10;
                                c11 += a10*b01 + a11*b11;
                        }

                        C[t] = c00;
                        C[t+1] = c01;
                        C[tt] = c10;
                        C[tt+1] = c11;
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
	int i , j , k;
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			register double sum=C[i*n+j];
			for(k=0;k<n;k++)
			{
				sum+=A[i*n+k]*B[k*n+j];
		
			}
			C[i*n+j]=sum;
		}
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i , j , k, iB , jB , kB;
	for(i=0;i<n;i+=BLOCK)
		for(j=0;j<n;j+=BLOCK)	
			for(k=0;k<n;k+=BLOCK)
			{
				for(iB=i;iB<i+BLOCK && iB<n;iB++)
					for(jB=j;jB<j+BLOCK && jB<n;jB++)
					{
						register double sum=C[iB*n+jB];
						for(kB=k;kB<k+BLOCK && kB<n;kB++)
							sum+=A[iB*n+kB]*B[kB*n+jB];
						C[iB*n+jB]=sum;
					}
			}
}

void jik(const double *A, const double *B, double *C, const int n) 
{
	int i , j , k;
	for(j=0;j<n;j++)
		for(i=0;i<n;i++)
		{
			register double sum=C[i*n+j];
			for(k=0;k<n;k++)
			{
				sum+=A[i*n+k]*B[k*n+j];
		
			}
			C[i*n+j]=sum;
		}

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i , j , k, iB , jB , kB;
	for(j=0;j<n;j+=BLOCK)
		for(i=0;i<n;i+=BLOCK)	
			for(k=0;k<n;k+=BLOCK)
			{
				for(jB=j;jB<j+BLOCK && jB<n;jB++)
					for(iB=i;iB<i+BLOCK && iB<n;iB++)
					{
						register double sum=C[iB*n+jB];
						for(kB=k;kB<k+BLOCK && kB<n;kB++)
							sum+=A[iB*n+kB]*B[kB*n+jB];
						C[iB*n+jB]=sum;
					}
			}
}

void kij(const double *A, const double *B, double *C, const int n) 
{
	int i , j , k;	
	for(k=0;k<n;k++)
		for(i=0;i<n;i++)
		{
			register double sum=A[i*n+k];
			for(j=0;j<n;j++)
			{
				C[i*n+j]+=sum*B[k*n+j];
		
			}
		}
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i , j , k, iB , jB , kB;
	for(k=0;k<n;k+=BLOCK)
		for(i=0;i<n;i+=BLOCK)	
			for( j=0;j<n;j+=BLOCK)
			{
				for(kB=k;kB<k+BLOCK && kB<n;kB++)
					for(iB=i;iB<i+BLOCK && iB<n;iB++)
					{
						register double sum=A[iB*n+kB];
						for(jB=j;jB<j+BLOCK && jB<n;jB++)
							C[iB*n+jB]+=sum*B[kB*n+jB];
					}
			}

}


void ikj(const double *A, const double *B, double *C, const int n) 
{
	int i , j , k;
	for(i=0;i<n;i++)
		for(k=0;k<n;k++)
		{
			register double sum=A[i*n+k];
			for(j=0;j<n;j++)
			{
				C[i*n+j]+=sum*B[k*n+j];
		
			}
		}
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i , j , k, iB , jB , kB;
	for(i=0;i<n;i+=BLOCK)
		for(k=0;k<n;k+=BLOCK)	
			for( j=0;j<n;j+=BLOCK)
			{
				for(iB=i;iB<i+BLOCK && iB<n;iB++)
					for(kB=k;kB<k+BLOCK && kB<n;kB++)
					{
						register double sum=A[iB*n+kB];
						for(jB=j;jB<j+BLOCK && jB<n;jB++)
							C[iB*n+jB]+=sum*B[kB*n+jB];
					}
			}


}

void jki(const double *A, const double *B, double *C, const int n) 
{
	int i , j , k;
	for(j=0;j<n;j++)
		for(k=0;k<n;k++)
		{
			register double sum=B[k*n+j];
			for(i=0;i<n;i++)
			{
				C[i*n+j]+=sum*A[i*n+k];
		
			}
		}
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i , j , k, iB , jB , kB;
	for(j=0;j<n;j+=BLOCK)
		for(k=0;k<n;k+=BLOCK)	
			for(i=0;i<n;i+=BLOCK)
			{
				for(jB=j;jB<j+BLOCK && jB<n;jB++)
					for(kB=k;kB<k+BLOCK && kB<n;kB++)
					{
						register double sum=B[kB*n+jB];
						for(iB=i;iB<i+BLOCK && iB<n;iB++)
							C[iB*n+jB]+=sum*A[iB*n+kB];
					}
			}
		
}

void kji(const double *A, const double *B, double *C, const int n) 
{
	int i , j , k;
	for(k=0;k<n;k++)
		for(j=0;j<n;j++)
		{
			register double sum=B[k*n+j];
			for(i=0;i<n;i++)
			{
				C[i*n+j]+=sum*A[i*n+k];
		
			}
		}
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i , j , k, iB , jB , kB;
	for(k=0;k<n;k+=BLOCK)
		for(j=0;j<n;j+=BLOCK)	
			for(i=0;i<n;i+=BLOCK)
			{
				for(kB=k;kB<k+BLOCK && kB<n;kB++)
					for(jB=j;jB<j+BLOCK && jB<n;jB++)
					{
						register double sum=B[kB*n+jB];
						for(iB=i;iB<i+BLOCK && iB<n;iB++)
							C[iB*n+jB]+=sum*A[iB*n+kB];
					}
			}
			
}

void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
	int i , j , k , iB, jB , kB;
        for( k=0;k<n;k+=b)
                for( i=0;i<n;i+=b)
                        for( j=0;j<n;j+=b)
                        {
                                for( kB=k;kB<k+b && kB<n;kB+=2)
                                        for( iB=i;iB<i+b && iB<n;iB+=2)
                                        {
                                                register int regA00=iB*n+kB;
                                                register int regA10=regA00+n;
                                                register double a00=A[regA00],a01=A[regA00+1],a10=A[regA10],a11=A[regA10+1];
                                                for(jB=j;jB<j+b && jB<n;jB+=2)
                                                {
                                                        register int regB00=kB*n+jB,regC00=iB*n+jB;
                                                        register int regB10=regB00+n,regC10=regC00+n;
                                                        register double b00=B[regB00],b01=B[regB00+1],b10=B[regB10],b11=B[regB10+1];
                                                        register double c00=C[regC00],c01=C[regC00+1],c10=C[regC10],c11=C[regC10+1];
                                                        C[regC00]+=a00*b00+a01*b10;
                                                        C[regC00+1]+=a00*b01+a01*b11;
                                                        C[regC10]+=a10*b00+a11*b10;
                                                        C[regC10+1]+=a10*b01+a11*b11;

                                                }
                                        }
                        }

}
