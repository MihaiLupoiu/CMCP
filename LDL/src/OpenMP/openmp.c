/*
 ============================================================================
 Name        : Assignment CMCP
 Author      : Mihaita Alexandru Lupoiu
 Version     : 0.0.1
 Description : LDL' Factorization using OpenMP and MPI
 ============================================================================
 */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "../tools.h"

#include <omp.h>

void ldl_Overwrite(double** A, size_t lA){
    int j,i,k;
    double* v = (double*) calloc(lA, sizeof(double));

    for (j = 0; j < lA ; ++j) {

        #pragma omp parallel for schedule(static) private(i)
        for (i = 0; i < j ; ++i) {
            v[i] = A[i][i]*A[j][i]; //Dependencia de datos, necesita tener la D (A[i][i]) y L (A[j][i]) anterior.
            //V es privada.
        }

        double ts = 0;
        #pragma omp parallel for schedule(static) reduction(+:ts)
        for (k = 0; k < j ; ++k) {
            ts = ts + A[j][k]*v[k]; // //Dependencia de datos, necesita tener la L (A[j][k]) anterior.
        }

        v[j] = A[j][j]-ts;
        A[j][j] = v[j]; // Se escribe la nueva D (A[j][j])

        #pragma omp parallel for schedule(static) private(k,ts,i)
        for (i = j+1; i < lA ; ++i) {
            ts = 0;
            for (k = 0; k < j ; ++k) {
                ts = ts + A[i][k]*v[k]; //Dependencia de datos, necesita tener la L (A[i][k]) anterior.
            }
            A[i][j] = (A[i][j]-ts)/v[j]; // // Se escribe la nueva L (A[i][j])
        }
    }
    free(v);
}

double frobeniusNorm(double **A, double **LDLT, size_t size1, size_t size2)
{
    double result = 0.0;
    int i,j;
    for(i = 0; i < size1; ++i)
    {
        for(j = 0; j < size2; ++j)
        {
            double value = A[i][j] - LDLT[i][j];
            result += value * value;
        }
    }
    return sqrt(result);
}

double norm(double** LDL, double** A, size_t size){

    double ** LD = make2dmatrix(size);

    int i,j,k = 0;

    for (i = 0; i < size; ++i) {
        for (j = 0; j <= i; j++) {
            if (i == j) {
                LD[i][j] = LDL[i][j];
            } else {
                LD[i][j] += LDL[i][j] * LDL[j][j];
            }
        }
    }

    double ** LT = make2dmatrix(size);
    for (i = 0; i < size; ++i) {
        for (j = i; j < size; j++) {
            if (i == j) {
                LT[i][j] = 1;
            } else {
                LT[i][j] = LDL[j][i];
            }

        }
    }

    double ** LDLT = make2dmatrix(size);
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; j++) {
            for (k = 0; k < size; k++) {
                LDLT[i][j] = LDLT[i][j] + LD[i][k]*LT[k][j];
            }
        }
    }

    double resultado = frobeniusNorm(A,LDLT, size,size);

    free2dmatrix(LD, size);
    free2dmatrix(LT, size);
    free2dmatrix(LDLT, size);

    return resultado;
}

int main(int argc, char *argv[]){

    int matrix_size,threads;
    int maxThreads = omp_get_max_threads();

    if(argc !=3){
        printf("Enter the size of matrix (N x N) where N = ");
        scanf("%d",&matrix_size);
        printf("Enter the CPUs number = ");
        scanf("%d",&threads);
    }
    else{
        matrix_size=atoi(argv[1]);
        threads=atoi(argv[2]);
    }

    if (threads >= maxThreads){
        threads = maxThreads;
    }
    omp_set_num_threads(threads);


    double **matrix=make2dmatrix(matrix_size);
    double **originalMatrix=make2dmatrix(matrix_size);

    initialize(matrix, matrix_size);
    copyMatrix(matrix, originalMatrix, matrix_size);

    //printMatrix(matrix,matrix_size);
    //printMatrix(originalMatrix,matrix_size);

/**
 * Code to Time the LDL' decompose
 */

    struct timeval t1, t2;
    double elapsedTime;

    //start timer
    gettimeofday(&t1, NULL);

    ldl_Overwrite(matrix, matrix_size);

    // stop timer
    gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec)+((t2.tv_usec - t1.tv_usec)/(1000.0*1000.0));

    //printMatrix(matrix,matrix_size);

    //printf("Frobenius Norm: %20.20f\n", norm(matrix, originalMatrix, matrix_size));

/*
    printf("\n**********************************\n\n");
    printf("Selected :%s\n","OpenMP");
    printf("Size of Matrix :%d \n",matrix_size);
    printf("DECOMPOSE TIME TAKEN : %f seconds\n",elapsedTime);
    printf("\n**********************************\n\n");
*/
    printf("%d;",matrix_size);
    printf("%f;",elapsedTime);
    printf("%d\n",threads);

    free2dmatrix(matrix,matrix_size);
    free2dmatrix(originalMatrix,matrix_size);
    return 0;
}
