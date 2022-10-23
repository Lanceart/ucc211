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
        
        if (max == 0) //expect condition
        {
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
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b) {
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    /* A, B and C are n x n matrices.
    /* This function computes C[:i,:j]+=A[:i,:k]*B[:k,:j] (the first i rows and k columuns of A multiplies the first k rows and j columuns of B added to the the first i rows and j columuns of C)
    /* b is the "block size" used in the dgemm.
    /* In fact this function won't be directly called in the tester code, so you can modify the declaration (parameter list) of mydgemm() if needed. 
    /* you may copy the code from the optimal() function or any of the other functions in your lab1 code (optimized code recommended).*/
    /* add your code here */
    
    int i1, j1, k1;

    for (i1 = i; i1 < i+b; i1 += 2) {
        for (k1 = k; k1 < k+b; k1 += 2) {
            // reg for A
            register int ra1 = i1*n + k1;// A[i1, k1]
            register int ra2 = ra1 + n;// A[(i1+1), k1]
            register double a00 = A[ra1];// A[i1, k1]
            register double a01 = A[ra1+1];// A[i1, k1+1]
            register double a10 = A[ra2];// A[(i1+1), k1]
            register double a11 = A[ra2+1];// A[(i1+1), (k1+1)]
            for (j1 = j; j1 < j+b; j1 += 2) {
                // reg for C
                register int rc1 = i1*n + j1;// C[i1, j1]
                register int rc2 = rc1 + n;// C[i1+1, j1]
                register double c00 = C[rc1];// C[i1, j1]
                register double c01 = C[rc1+1];// C[i1, j1+1]
                register double c10 = C[rc2];// C[i1+1, j1]
                register double c11 = C[rc2+1];// C[i1+1, j1+1]
                // reg for B
                register int rb1 = k1*n + j1;// B[k1, j1]
                register int rb2 = rb1 + n;// B[(k1+1), j1]
                register double b00 = B[rb1];// B[k1, j1]
                register double b01 = B[rb1+1];// B[k1, j1+1]
                register double b10 = B[rb2];// B[(k1+1), j1]
                register double b11 = B[rb2+1];// B[k1+1, j1+1]
                // calculate C
                c00 -= a00 * b00 + a01 * b10;
                c01 -= a00 * b01 + a01 * b11;
                c10 -= a10 * b00 + a11 * b10;
                c11 -= a10 * b01 + a11 * b11;
                C[rc1] = c00;
                C[rc1+1] = c01;
                C[rc2] = c10;
                C[rc2+1] = c11;
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
int mydgetrf_block(double *A, int *ipiv, int n, int b) {
    int i, j, k;
    int max_Aii;
    int ib;
    double max;
    double* temp_Aii = (double*)malloc(n * sizeof(double)); // a temp row for swap
    
    for (ib = 0;ib < n;ib += b) {
        for (i = ib; i < ib+b; i++) {
            max = fabs(A[i*n + i]);
            max_Aii = i;
            for (j = i+1;j < n;j++) {
                if ( fabs(A[j*n + i]) > max ) {
                    max = fabs(A[j*n + i]);
                    max_Aii = j; // |A[j,i]| is the biggest in rest of i-th col
                }
            }
            // swap row i and max_Aii and the corrosponding vector
            if (max == 0) {
                return -1;
            } else {
                if (max_Aii != i) { 
                    // swap i-th and max_Aii-th element of vector pivot
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[max_Aii];
                    ipiv[max_Aii] = temp;
                    // swap i-th row and max_Aii-th row
                    memcpy(temp_Aii, A + i*n, n*sizeof(double));
                    memcpy(A + i*n, A + max_Aii*n, n*sizeof(double));
                    memcpy(A + max_Aii*n, temp_Aii, n*sizeof(double));
                }
            }
            for (j = i+1; j < n; j++) {
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                for (k = i+1; k < ib+b; k++) {
                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }
            }
        }
        // A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n), update next b rows of U
        for (i = ib; i < ib+b; i++) {
            for (j = ib+b;j < n;j++) { // j:[end+1=ib+b-1+1, n]
                double sum = 0.0;
                for (k = ib;k < i;k++) {
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }

        // update A(end+1=ib+b-1+1 : n, end+1:n) block by block
        for (i = ib+b;i < n;i+=b) {
            for (j = ib+b;j < n;j+=b) {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }

    return 0;
}
