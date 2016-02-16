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
#include <math.h>
#include "../tools.h"
#include <sys/time.h>



void print_secuencial_ldl(double** A, size_t lA ){

    int j,i,k;
    double* v = (double*) calloc(lA, sizeof(double));

    for (j = 0; j < lA ; ++j) {
        printf("j=%d\n",(int)j);

        for (i = 0; i <= j-1 ; ++i) {
            printf("\ti=%d <= %d\n",(int)i,(int)j-1);
            v[i] = A[i][i]*A[j][i];
            printf("\t\tv[%d]=D[%d]*L[%d][%d]\n",(int)i,(int)i,(int)j,(int)i);
        }

        double ts = 0;
        printf("\tts=%d\n",(int)ts);
        for (k = 0; k <= j-1 ; ++k) {
            printf("\tk=%d <= %d\n",(int)k,(int)j-1);
            ts = ts + A[j][k]*v[k];
            printf("\t\tts=ts+L[%d][%d]*v[%d]\n",(int)j,(int)k,(int)k);
        }

        v[j] = A[j][j]-ts;
        printf("\tv[%d]=A[%d][%d]-ts\n",(int)j,(int)j,(int)j);
        printf("\tD[%d]=v[%d]\n",(int)j,(int)j);

        for (i = j+1; i < lA ; ++i) {
            printf("\ti=%d < %d \n",(int)i,(int)lA);
            ts = 0;
            printf("\t\tts=%d\n",(int)ts);
            for (k = 0; k <= j-1 ; ++k) {
                printf("\t\tk=%d < %d\n",(int)k,(int)j-1);
                ts = ts + A[i][k]*v[k];
                printf("\t\t\tts=ts+L[%d][%d]*v[%d]\n",(int)i,(int)k,(int)k);
            }
            A[i][j] = (A[i][j]-ts)/v[j];
            printf("\t\tL[%d][%d]=(A[%d][%d]-ts)/v[%d]\n",(int)i,(int)j,(int)i,(int)j,(int)j);

        }
    }
}

void ldl_Overwrite(double** A, size_t lA){
    int j,i,k;
    double* v = (double*) calloc(lA, sizeof(double));

    for (j = 0; j < lA ; ++j) {

        for (i = 0; i < j ; ++i) {
            v[i] = A[i][i]*A[j][i]; //Dependencia de datos, necesita tener la D (A[i][i]) y L (A[j][i]) anterior.
                                    //V es privada.
        }

        double ts = 0;
        for (k = 0; k < j ; ++k) {
            ts = ts + A[j][k]*v[k]; // //Dependencia de datos, necesita tener la L (A[j][k]) anterior.
        }

        v[j] = A[j][j]-ts;
        A[j][j] = v[j]; // Se escribe la nueva D (A[j][j])

        for (i = j+1; i < lA ; ++i) {
            ts = 0;
            for (k = 0; k < j ; ++k) {
                ts = ts + A[i][k]*v[k]; //Dependencia de datos, necesita tener la L (A[i][k]) anterior.
            }
            A[i][j] = (A[i][j]-ts)/v[j]; // // Se escribe la nueva L (A[i][j])
        }
     }
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
            //printf("%20.18f - %20.18f\n",A[i][j], LDLT[i][j]);

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
/*
    printf("=============================================================\n");
    printf("                           LD\n");
    printf("=============================================================\n");
    printMatrix(LD, size);
*/
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

    /*
    printf("=============================================================\n");
    printf("                           LT\n");
    printf("=============================================================\n");
    printMatrix(LT, size);
    */
    double ** LDLT = make2dmatrix(size);
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; j++) {
            for (k = 0; k < size; k++) {
                LDLT[i][j] = LDLT[i][j] + LD[i][k]*LT[k][j];
            }
        }
    }
/*
    printf("=============================================================\n");
    printf("                           LDLT\n");
    printf("=============================================================\n");
    printMatrix(LDLT, size);
*/
    double resultado = frobeniusNorm(A,LDLT, size,size);

    free2dmatrix(LD, size);
    free2dmatrix(LT, size);
    free2dmatrix(LDLT, size);

    return resultado;
}

int main(int argc, char *argv[]){


    int matrix_size;

    if(argc !=2){
        printf("Enter the size of matrix (N x N) where N = ");
        scanf("%d",&matrix_size);
    }
    else{
        matrix_size=atol(argv[1]);
    }

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
    //print_secuencial_ldl(matrix, matrix_size);

    // stop timer
    gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec)+((t2.tv_usec - t1.tv_usec)/(1000.0*1000.0));

    //printMatrix(matrix,matrix_size);

    //printf("Frobenius Norm: %20.20f\n", norm(matrix, originalMatrix, matrix_size));

/*
    printf("\n**********************************\n\n");
    printf("Selected :%s\n","Sequential");
    printf("Size of Matrix :%d \n",matrix_size);
    printf("DECOMPOSE TIME TAKEN : %f seconds\n",elapsedTime);
    printf("\n**********************************\n\n");
*/

    printf("%d;",matrix_size);
    printf("%f\n",elapsedTime);

    free2dmatrix(matrix,matrix_size);
    return 0;
}
