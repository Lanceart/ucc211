#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
	int i,j,k;
	for (i=0;i<n;i++)
	{
		for(j=0;j<n;j++){
			for(k=0;k<n;k++)
			C[n*i+j]+=A[i*n+k]*B[k*n+j];
			}
	}
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
	int i,k,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
		{
			register double temps = C[i*n+j];
			for(k=0;k<n;k++){
				temps+=A[i*n+k]*B[k*n+j];
				}
			C[i*n+j] = temps;
			}
	}
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
int i, j, k;
	for (i = 0; i < n; i+=2){
		for (j = 0; j < n; j+=2){
			for (k = 0; k < n; k+=2){
				C[i*n + j]= A[i*n + k]*B[k*n + j] + A[i*n + k+1]*B[(k+1)*n + j]+ C[i*n + j];
				C[(i+1)*n + j]= A[(i+1)*n + k]*B[k*n + j] + A[(i+1)*n +
k+1]*B[(k+1)*n + j]+ C[(i+1)*n + j];
				C[i*n + (j+1)]= A[i*n + k]*B[k*n + (j+1)] + A[i*n +
k+1]*B[(k+1)*n + (j+1)]+ C[i*n + (j+1)];
				C[(i+1)*n + (j+1)]= A[(i+1)*n + k]*B[k*n + (j+1)]+ A[(i+1)*n + k+1]*B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
	}
	}
	}
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 2)
    {
        for (j = 0; j < n; j +=2 )
        {
            register double C_00 = C[i * n + j];
            register double C_01 = C[i * n + (j + 1)];
            register double C_10 = C[(i + 1) * n + j];
            register double C_11 = C[(i + 1) * n + (j + 1)];

            for (k = 0; k < n; k+=2)
            {
		register int a_1 = i*n+k;
		register int a_2 = a_1+n;
		register int b_1 = k*n+j;
		register int b_2 = b_1+n;

                register double A_00 = A[a_1];
                register double A_01 = A[a_1+1];
		register double A_10 = A[a_2];
		register double A_11 = A[a_2+1];

                register double B_00 = B[b_1];
                register double B_01 = B[b_1+1];
		register double B_10 = B[b_2];
		register double B_11 = B[b_2+1];
		
		C_00 += A_00*B_00 + A_01*B_10;
		C_01 += A_00*B_01 + A_01*B_11;
		C_10 += A_10*B_00 + A_11*B_10;
		C_11 += A_10*B_01 + A_11*B_11;		
            }
         	
		C[i*n+j]= C_00;C[i*n+j+1]=C_01;
		C[(i+1)*n+j]=C_10;C[(i+1)*n+j+1]=C_11;	   
        }
    }
}

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            register double R = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                R += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = R;
        }
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    int i, k, j, i1, k1, j1;
    for (i = 0; i < n; i += b)
    {
        for (j = 0; j < n; j += b)
        {
            for (k = 0; k < n; k += b)
            {
                for (i1 = i; i1 < i + b && i1 < n; i1++)
                {
                    for (j1 = j; j1 < j + b && j1 < n; j1++)
                    {
                        register double R = C[i1 * n + j1];
                        for (k1 = k; k1 < k + b && k1 < n; k1++)
                        {
                            R += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = R;
                    }
                }
            }
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n) 
{
	    int j, i, k;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            register double R = C[i * n + j];
            for (k = 0; k < n; k++)
            {
                R += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = R;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{

    int j, i, k, j1, i1, k1;
    for (j = 0; j < n; j += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (k = 0; k < n; k += b)
            {
                for (j1 = j; j1 < j + b && j1 < n; j1++)
                {
                    for (i1 = i; i1 < i + b && i1 < n; i1++)
                    {
                        register double R = C[i1 * n + j1];
                        for (k1 = k; k1 < k + b && k1 < n; k1++)
                        {
                            R += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = R;
                    }
                }
            }
        }
    }
}


void kij(const double *A, const double *B, double *C, const int n) 
{
	    int k, i, j;
    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            register double R = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += R * B[k * n + j];
            }
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
	    int k, i, j, k1, i1, j1;
    for (k = 0; k < n; k += b)
    {
        for (i = 0; i < n; i += b)
        {
            for (j = 0; j < n; j += b)
            {
                for (k1 = k; k1 < k + b && k1 < n; k1++)
                {
                    for (i1 = i; i1 < i + b && i1 < n; i1++)
                    {
                        register double R = A[i1 * n + k1];
                        for (j1 = j; j1 < j + b && j1 < n; j1++)
                        {
                            C[i1 * n + j1] += R * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
	int i, k, j;
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < n; k++)
        {
            register double R = A[i * n + k];
            for (j = 0; j < n; j++)
            {
                C[i * n + j] += R * B[k * n + j];
            }
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{

int i, k, j, i1, k1, j1;
    for (i = 0; i < n; i += b)
    {
        for (k = 0; k < n; k += b)
        {
            for (j = 0; j < n; j += b)
            {
                for (i1 = i; i1 < i + b && i1 < n; i1++)
                {
                    for (k1 = k; k1 < k + b && k1 < n; k1++)
                    {
                        register double R = A[i1 * n + k1];
                        for (j1 = j; j1 < j + b && j1 < n; j1++)
                        {
                            C[i1 * n + j1] += R * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }
}

void jki(const double *A, const double *B, double *C, const int n) 
{
int j, k, i;
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k++)
        {
            register double R = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * R;
            }
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{

int j, k, i, j1, k1, i1;
    for (j = 0; j < n; j += b)
    {
        for (k = 0; k < n; k += b)
        {
            for (i = 0; i < n; i += b)
            {
                for (j1 = j; j1 < j + b && j1 < n; j1++)
                {
                    for (k1 = k; k1 < k + b && k1 < n; k1++)
                    {
                        register double R = B[k1 * n + j1];
                        for (i1 = i; i1 < i + b && i1 < n; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * R;
                        }
                    }
                }
            }
        }
    }
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    int k, j ,i;
    for (k = 0; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            register double R = B[k * n + j];
            for (i = 0; i < n; i++)
            {
                C[i * n + j] += A[i * n + k] * R;
            }
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{

    int k, j, i, k1, j1, i1;
    for (k = 0; k < n; k += b)
    {
        for (j = 0; j < n; j += b)
        {
            for (i = 0; i < n; i += b)
            {
                for (k1 = k; k1 < k + b && k1 < n; k1++)
                {
                    for (j1 = j; j1 < j + b && j1 < n; j1++)
                    {
                        register double R = B[k1 * n + j1];
                        for (i1 = i; i1 < i + b && i1 < n; i1++)
                        {
                            C[i1 * n + j1] += A[i1 * n + k1] * R;
                        }
                    }
                }
            }
        }
    }
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}
