#include "../include/for_you_to_do.h"

#include <math.h>

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/



int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, j, k, max_index, tmp2;
    double max, tmp1;
    // used for swapping rows
    double * tmp_row = (double *)malloc(sizeof(double) * n);
    
    for (i = 0; i < n; i++)
    {
        max_index = i;
        max = fabs(A[i * n + i]);
        for (j = i + 1; j < n; j++)
        {
            tmp1 = fabs(A[j * n + i]);
            if (max < tmp1)
            {
                max_index = j;
                max = tmp1;
            }
        }
        
        // if the matrix is singular
        if (max == 0)
        {
	    perror("LU factorization failed: coefficient matrix is singular.\n");
            return -1;
        }
           
        if (max_index != i)
        {	
            tmp2 = ipiv[i];
            ipiv[i] = ipiv[max_index];
            ipiv[max_index] = tmp2;
		
            // swap rows of A
            memcpy(tmp_row, A + i * n, sizeof(double) * n);
            memcpy(A + i * n, A + max_index * n, sizeof(double) * n);
            memcpy(A + max_index * n, tmp_row, sizeof(double) * n);
        }
    
        for (j = i + 1; j < n; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++)
            {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
    free(tmp_row);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i, j;
    double sum = 0;
    double * tmp_B = (double *)malloc(sizeof(double) * n);
    
    if (UPLO == 'L')
    {
	   
        for (i = 0; i < n; i++)
        {
            tmp_B[i] = B[ipiv[i]];
        }
        
        for (i = 0; i < n; i++)
        {
            sum = tmp_B[i];
            for (j = 0; j < i; j++)
            {
                sum -= B[j] * A[i * n + j];
            }
            B[i] = sum;
        }
    }
    else if (UPLO == 'U')
    {
	// Ux = y, backward subtitution     
        for (i = n - 1; i >= 0; i--)
        {
            sum = 0;
            for (j = i + 1; j < n; j++)
            {
                sum += B[j] * A[i * n + j];
            }
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    free(tmp_B);
    return;
}


int get_block_size(){
    //return the block size you use in your matrix multiplication code.
    /*add your code here, 128 is an example and can be modified. */
    // With the returned block size the test code will adaptively adjust the input sizes to avoid corner cases.
    return 128;
  
}

//The matrix multiplication function used in blocked GEPP.
// You need to let the mydgemm adapt to non-square inputs.
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* A, B and C are n x n matrices.
    /* This function computes C[:i,:j]+=A[:i,:k]*B[:k,:j] (the first i rows and k columuns of A multiplies the first k rows and j columuns of B added to the the first i rows and j columuns of C)
    /* b is the "block size" used in the dgemm.
    /* In fact this function won't be directly called in the tester code, so you can modify the declaration (parameter list) of mydgemm() if needed. 
    /* you may copy the code from the optimal() function or any of the other functions in your lab1 code (optimized code recommended).*/
    /* add your code here */
    	int i , j , k , iB, jB , kB;
        for( k = 0;k < maty ; k += b)
                for( i = 0;i < matx;i += b)
                        for( j = 0;j < matx ;j += b)
                        {
                                for( kB = k;kB < k + b && kB < maty;kB += 2)
                                        for( iB = i;iB <i + b && iB < matx;iB += 2)
                                        {
                                                register int regA00 = iB *n + kB;
                                                register int regA10 = regA00 + n;
                                                register double a00 = A[regA00],a01 = A[regA00 + 1],a10 = A[regA10],a11 = A[regA10 + 1];
                                                for(jB = j;jB < j + b && jB < matx;jB += 2)
                                                {
                                                        register int regB00 = kB * n + jB,regC00 = iB * n + jB;
                                                        register int regB10 = regB00 + n,regC10 = regC00 + n;
                                                        register double b00 = B[regB00],b01 = B[regB00 + 1],b10 = B[regB10],b11 = B[regB10 + 1];
                                                        register double c00 = C[regC00],c01 = C[regC00 + 1],c10 = C[regC10],c11 = C[regC10 + 1];
                                                        C[regC00] -= a00 * b00 + a01 * b10;
                                                        C[regC00+1] -= a00 * b01 + a01 * b11;
                                                        C[regC10] -= a10 * b00 + a11 * b10;
                                                        C[regC10+1] -= a10 * b01 + a11 * b11;
                                                }
                                        }
                        }
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int i, maxind, k, j, ib, end, temps, t;
		double max;
		for(ib = 0; ib < n - 1; ib += b)
		{
			/*Partial Pivoting*/
			end = ((n-1) > (ib + b -1)) ? (ib + b - 1) : n-1;
			//printf("end = %d\n",end);
			for(i = ib; i <= end; i++)
			{
				maxind = i;
				max = fabs(A[i*n+i]);
				for(k = i+1; k < n; k++)
				{
					if(fabs(A[k * n + i]) > max)
					{
						maxind = k;
						max = fabs(A[k * n + i]);
					}
				}
				if(max == 0) return -1;
				else if (maxind !=i)
				{

					/*Save Pivot Infortmation*/
					temps = ipiv[i];
					ipiv[i] = ipiv[maxind];
					ipiv[maxind] = temps;
					/*Swap rows*/
					for(j = 0; j < n; j++)
					{
						double tempv;
						tempv = A[i * n + j];
						A[i * n + j] = A[maxind * n + j];
						A[maxind * n + j] = tempv;
					}
				}

				/*Update columns i+1 to end*/
				for(j = i + 1; j < n; j++)
				{
					A[j * n + i] = (double)A[j * n + i] / A[i * n + i];
					for(t = i + 1;t <= end; t++)
					{
						A[j*n+t] = A[j*n+t] - A[j*n+i] * A[i*n+t];
					}
				}
			}

			/*inv(LL)*/
			/*double y;y = (double *) malloc(sizeof(double) * (end - ib + 1) * (n - end));y[0] =; */
			for(i = ib; i <= end; i++)
			{
				for(k = end +1; k < n; k++)
				{
					double sum = 0;
					for(j = ib; j < i; j++)
					{
						sum += A[i * n + j] * A [j * n + k];
					}
					A[i * n + k] -= sum;
				}
			}
			/*Delayed update of rest of matrix using matrix-matrix multiplication*/
			/*void mydgemm(double *A, double *B, double *C, int n, int matx, int maty, int b)*/
			//if(end!=n)
			mydgemm(&A[(end+1) * n + ib], &A[ib * n + end +1], &A[(end+1) * n + (end + 1)], n , (n - end - 1) , (end-ib+1/*=b*/), 32);
		}	

    return 0;
}
