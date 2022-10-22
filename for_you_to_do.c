#include "../include/for_you_to_do.h"



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
    /* add your code here */
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
    return 0;
}

