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
    double * temps_double = (double *)malloc(sizeof(double) * n);
    
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
		
        //
        /*
        memcpy(temps_double, A + max_index * n, sizeof(double) * n);
            memcpy(A + max_index * n, A + A + i * n, sizeof(double) * n);
            memcpy(A + A + i * n, temps_double, sizeof(double) * n);
        */
            //A + max_index * n
            memcpy(temps_double, A + max_index * n, sizeof(double) * n);
            memcpy(A + max_index * n, A + A + i * n, sizeof(double) * n);
            memcpy(A + A + i * n, temps_double, sizeof(double) * n);
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
    free(temps_double);
    return 0;
}


void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i, j;
    double sum = 0;
    double * tmp_B = (double *)malloc(sizeof(double) * n);
    
    if (UPLO == 'U')
    {   
        for (i = n - 1; i >= 0; i--)
        {
            sum = 0;
            int temp_index = i*n;
            for (j = i + 1; j < n; j++)
            {
                sum += B[j] * A[temp_index + j];
            }
            B[i] = (B[i] - sum) / A[temp_index + i];
        }
    }


    if (UPLO == 'L')
    {
	   
        for (i = 0; i < n; i++)
        {
            tmp_B[i] = B[ipiv[i]];
        }
        
        for (i = 0; i < n; i++)
        {
            sum = tmp_B[i];
            int temp_index = i*n;
            for (j = 0; j < i; j++)
            {
                sum -= B[j] * A[temp_index + j];
            }
            B[i] = sum;
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
            register int t_temp = i1*n + k1;
            register int q_temp = t_temp + n;


            register double A_00 = A[t_temp];
            register double A_01 = A[t_temp+1];
            register double A_10 = A[q_temp];
            register double A_11 = A[q_temp+1];
            for (j1 = j; j1 < j+b; j1 += 2) {
                
                register int rb1 = k1*n + j1;
                register int rb2 = rb1 + n;
                register double B_00 = B[rb1];
                register double B_01 = B[rb1+1];
                register double B_10 = B[rb2];
                register double B_11 = B[rb2+1];



                register int rc1 = i1*n + j1;
                register int rc2 = rc1 + n;
                register double C_00 = C[rc1];
                register double C_01 = C[rc1+1];
                register double C_10 = C[rc2];
                register double C_11 = C[rc2+1];



                C[rc1] = C_00 - ( A_00 * B_00 + A_01 * B_10);
                C[rc1+1] = C_01 - (A_00 * B_01 + A_01 * B_11);
                C[rc2] = C_10 - (A_10 * B_00 + A_11 * B_10);
                C[rc2+1] = C_11 - (A_10 * B_01 + A_11 * B_11);
            }
        }
    }
    return;
}



int mydgetrf_block(double *A, int *ipiv, int n, int b) {
    double* temps_index_double = (double*)malloc(n * sizeof(double)); 
    int max_index_me,temp_i,i, j, k;
    double max;
    
    
    for (temp_i = 0;temp_i < n;temp_i += b) {
        //if(temp_i == n)break;
        for (i = temp_i; i < temp_i+b; i++) {
            max = fabs(A[i*n + i]);
            max_index_me = i;
            for (j = i+1;j < n;j++) {
                if ( fabs(A[j*n + i]) > max ) {
                    max = fabs(A[j*n + i]);
                    max_index_me = j; 
                }
            }

            if (max == 0) {
                return -1;
            } else {
                if (max_index_me != i) { 

                    int temp = ipiv[i];
                    ipiv[i] = ipiv[max_index_me];
                    ipiv[max_index_me] = temp;

                    memcpy(temps_index_double, A + i*n, n*sizeof(double));
                    memcpy(A + i*n, A + max_index_me*n, n*sizeof(double));
                    memcpy(A + max_index_me*n, temps_index_double, n*sizeof(double));

                    /*
                    memcpy(temps_index_double, A + i*n, n*sizeof(double));
                    memcpy(A + i*n, A + max_index_me*n, n*sizeof(double));
                    memcpy(A + max_index_me*n, temps_index_double, n*sizeof(double));
                    */
                }
            }

            if(max ==0){return -1;}
            for (j = i+1; j < n; j++) {
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                for (k = i+1; k < temp_i+b; k++) {
                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }
            }
        }
        double sum = 0.0;
        for (i = temp_i; i < temp_i+b; i++) {
            for (j = temp_i+b;j < n;j++) { 
                sum = 0.0;
                for (k = temp_i;k < i;k++) {
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }


        for (i = temp_i+b;i < n;i+=b) {
            for (j = temp_i+b;j < n;j+=b) {
                mydgemm(A, A, A, n, i, j, temp_i, b);
            }
        }
    }

    return 0;
}
