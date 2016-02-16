/*
 ============================================================================
 Name        : Assignment CMCP
 Author      : Mihaita Alexandru Lupoiu
 Version     : 0.0.1
 Description : LDL' Factorization using OpenMP and MPI
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include "tools.h"

void printMatrix(double** matrix,size_t size){
    printf("\n *************** MATRIX ****************\n\n");
    int i,j;
    for(i = 0; i < size; i++) {
        for (j = 0; j < size; ++j) {
            printf(" %f ",matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printMatrixCustom(double* matrix,size_t size_r,size_t size_c){
    printf("\n *************** MATRIX ****************\n\n");
    int i,j;
    for(i = 0; i < size_r; i++) {
        for (j = 0; j < size_c; ++j) {
            printf(" %f ",matrix[i*size_c+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printMatrixCustomCPU(double* matrix,size_t size_r,size_t size_c, int poces){
    printf("\n *************** MATRIX Process %d ****************\n\n",poces);
    int i,j;
    for(i = 0; i < size_r; i++) {
        for (j = 0; j < size_c; ++j) {
            printf(" %f ",matrix[i*size_c+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printVector(double* vector, size_t size){
    printf("\n *************** VECTOR ****************\n\n");
    int i;
    for(i = 0; i < size; i++) {
        printf(" %f ",vector[i]);
    }
    printf("\n");
}

void initialize(double **matrix ,size_t size){
    //double A[3][3] = {{2,2,1},{2,5,2},{1,2,2}};

    int i,j, randomValue;
    for(i = 0; i < size; i++) {
        for (j = 0; j <= i; j++) {
            if(i==j){
                //matrix[i][j] = A[i][j];
                matrix[i][j]= (rand() % 100 + 1)+200;
            }else{
                //matrix[i][j] = A[i][j];
                //matrix[j][i] = A[j][i];
                randomValue = rand() % 100 + 1;
                matrix[i][j] = randomValue;
                matrix[j][i] = randomValue;
            }
        }
    }
}

void initializeCust(double *matrix ,size_t size){
    //double A[3][3] = {{2,2,1},{2,5,2},{1,2,2}};

    int i,j, randomValue;
    for(i = 0; i < size; i++) {
        for (j = 0; j <= i; j++) {
            if(i==j){
                //matrix[i+j*size] = A[i][j];

                matrix[i+j*size] = (rand() % 100 + 1)+200;
            }else{
                //matrix[i+j*size] = A[i][j];
                //matrix[j+i*size] = A[j][i];

                randomValue = rand() % 100 + 1;
                matrix[i+j*size] = randomValue;
                matrix[j+i*size] = randomValue;
            }
        }
    }
}

void copyMatrix(double **matrix, double **matrixToCopy ,size_t size){
    int i,j;
    for(i = 0; i < size; i++) {
        for (j = 0; j <= i; j++) {
                matrixToCopy[i][j] = matrix[i][j];
                matrixToCopy[j][i] = matrix[j][i];
            }
    }
}
void copyMatrixCustom(double *matrix, double *matrixToCopy ,size_t size){
    int i,j;
    for(i = 0; i < size; i++) {
        for (j = 0; j <= i; j++) {
            matrixToCopy[i*size+j] = matrix[i*size+j];
            matrixToCopy[j*size+i] = matrix[j*size+i];
        }
    }
}

double **make2dmatrix(size_t size) {
    int i;
    double **m;
    m = (double**)calloc(size,sizeof(double*));
    for (i=0;i<size;i++)
        m[i] = (double*)calloc(size,sizeof(double));
    return m;
}

void free2dmatrix(double **matrix, size_t size) {
    int i;
    if (!matrix) return;
    for(i=0;i<size;i++)
        free(matrix[i]);
    free(matrix);
}

double **make2dmatrixCustom(size_t size_r,size_t size_c){
    int i;

    double *dado = (double *)calloc(size_r*size_c,sizeof(double));

    double **array = (double **)calloc(size_r,sizeof(double*));

    for(i = 0; i < size_r; i++)
        array[i] = &(dado[size_c*i]);

    return array;
};

/* muestra de tiempo
    struct timeval t1, t2;
    double elapsedTime;

    //start timer
    //gettimeofday(&t1, NULL);

    //Codigo

   // stop timer
    //gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    //elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
    //printf(" --WFS PASADO %f -- \n", elapsedTime);

 */



/*
int main(int argc, char *argv[]){
    int matrix_size;

    if(argc !=2){
        printf("Enter the size of matrix (N x N) where N = ");
        scanf("%lu",&matrix_size);
    }
    else{
        matrix_size=atol(argv[1]);
    }

    double **matrix=make2dmatrix(matrix_size);

    initialize(matrix, matrix_size);

    printmatrix(matrix,matrix_size);
*/
    /**
     * Code to Time the LDL' decompose
     */
/*    clock_t startTimer, stopTimer;
    double time_spent;
    startTimer = clock();

    //decomposeSerial(matrix,matrix_size);

    stopTimer = clock();
    time_spent = ((double)(stopTimer - startTimer)) / CLOCKS_PER_SEC;

    printmatrix(matrix,matrix_size);

    printf("\n**********************************\n\n");
    printf("Selected :%s\n","Sequential");
    printf("Size of Matrix :%lu \n",matrix_size);
    printf("DECOMPOSE TIME TAKEN : %f seconds\n",time_spent);
    printf("\n**********************************\n\n");

    free2dmatrix(matrix,matrix_size);
    return 0;
}
*/