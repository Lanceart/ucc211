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
    register int i1, j1, k1; 
    // used for optimizing	  
    register int in1, in2, in3, j2, j3, kn;	
    register int n1 = i + b > n ? n : i + b;
    register int n2 = j + b > n ? n : j + b;
    register int n3 = k + b > n ? n : k + b;

    for (i1 = i; i1 < n1; i1 += 3)
    {
        for (j1 = j; j1 < n2; j1 += 3)
        {
	    in1 = i1 * n;
	    in2 = (i1 + 1) * n;
	    in3 = (i1 + 2) * n;
	    j2 = j1 + 1;
	    j3 = j1 + 2;	
            register double C_0_0 = C[in1 + j1];
            register double C_0_1 = C[in1 + j2];
            register double C_0_2 = C[in1 + j3];
            register double C_1_0 = C[in2 + j1];
            register double C_1_1 = C[in2 + j2];
            register double C_1_2 = C[in2 + j3];
            register double C_2_0 = C[in3 + j1];
            register double C_2_1 = C[in3 + j2];
            register double C_2_2 = C[in3 + j3];
		
	    for (k1 = k; k1 < n3; k1++)
            {
		kn = k1 * n + j1;    
                register double A_0 = A[in1 + k1];
                register double A_1 = A[in2 + k1];
                register double A_2 = A[in3 + k1];
                register double B_0 = B[kn];
                register double B_1 = B[kn + 1];
                register double B_2 = B[kn + 2];

                C_0_0 -= A_0 * B_0;
                C_0_1 -= A_0 * B_1;
                C_0_2 -= A_0 * B_2;
                C_1_0 -= A_1 * B_0;
                C_1_1 -= A_1 * B_1;
                C_1_2 -= A_1 * B_2;
                C_2_0 -= A_2 * B_0;
                C_2_1 -= A_2 * B_1;
                C_2_2 -= A_2 * B_2;
            }
		
            C[in1 + j1] = C_0_0;
            C[in1 + j2] = C_0_1;
            C[in1 + j3] = C_0_2;
            C[in2 + j1] = C_1_0;
            C[in2 + j2] = C_1_1;
            C[in2 + j3] = C_1_2;    
            C[in3 + j1] = C_2_0;
            C[in3 + j2] = C_2_1;
            C[in3 + j3] = C_2_2;
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
    int i,j,k,ic,t, maxind;
    double max;

    for(ic = 0; ic <n - 1;ic +=b){
        // if(ic > n){
        //     printf("error in ic\n\n\n");
        // }

        for(i = ic; i < ic+b ; i++){
            // if(i > n){
            //     printf("error in i\n\n\n");
            // }

            // pivoting the matrix to find the maximum number
            maxind = i;
            max = fabs(A[i*n + i]);
            for(t = i+1; t < n; t++){

                if(fabs(A[t*n + i]) > max){
                    maxind = t;
                    max = fabs(A[t*n + i]);
                }
            }

            
            if(max == 0){
                printf("LU factoration failed: coefficient matrix is singular");
                return -1;
            }
            else{

                //The case that need to swap the row
                if(maxind != i){
                    //swap pivoting information
                    int temps= ipiv[i];
                    ipiv[i] = ipiv[maxind];
                    ipiv[maxind] = temps;

                    //swap row for matrix method 1
                    // int j;
                    // for(j = 0; j < n; j++){
                    //     double k;
                    //     k = A[i * n + j];
                    //     A[i * n + j] = A[maxind * n + j];
                    //     A[maxind * n + j] = k;
                    // }

                    //swap row method 2 -- seem like not too much difference than method one, but this look a bit cleaner
                    double trow[n];
                    memcpy(trow, A + i * n, n*sizeof(double));
                    memcpy(A + i * n, A + maxind * n, n*sizeof(double));
                    memcpy(A + maxind * n, trow, n*sizeof(double));
                }

            }

            //update the lower triangle of A(ic:end , ic:end) and A(end+1:n , ic:end)
            for(j = i + 1; j <n;j++){
                // if(j > n){
                //     printf("error in j\n\n\n");
                // }

                A[j*n + i] = A[j*n + i] / A[i*n + i];

                //block version
                for(k = i + 1; k < ic + b; k++){

                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }

                //naive version - to test the top part of the code work
                // for(k = i + 1; k < n; k++){
                //     A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
                // }
            }
        }


        //update A(ic:end, end+1:n), basically same method as before, use the value store in A(ic:n, ic:end)
        register double total;
        //end = ic + b
        for(i = ic; i < ic + b; i++){
            // if(i > n){
            //     printf("error in i\n\n\n");
            // }

            for(j= ic + b;j < n;j++){
                // if(j > n){
                //     printf("error in j\n\n\n");
                // }

                total = 0;
                for(k = ic; k < i; k++){
                    // if(k > n){
                    //     printf("error in k\n\n\n");
                    // }

                    //naive version, abandon
                    // A[i*n - j] -= A[i*n + k] * A[k*n + j];

                    //new version, reduce access element
                    total += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= total;
            }
        }

        // update A(end + 1: n , end + 1 : n)
        // end = ic + b
        mydgemm(A, A, A,n, ic + b, ic + b, ic, b);
    }
return 0;
}
