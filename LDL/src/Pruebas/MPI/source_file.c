//
// Created by Mihaita Alexandru Lupoiu on 10/11/15.
//

/* C Example */
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <malloc.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <time.h>
//#include <omp.h>

void ctimer(double*, double*, double*);

void iniMatrix(double *A, int n) {
    int j, i;
    for (j=0; j<n; j++){
        for (i=0; i<n; i++){
            A[i+j*n] = ((double) rand()/RAND_MAX);
        }
    }
}

void printC(double *A, int n, int m){
    int i,j;
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++)
            printf("%5.2f ",A[i+j*n]);
        printf("\n");
    }
}

void vectMatrix(double *A, double *B, double *C, int n, int spl) {
    int i, j, k;
#pragma omp parallel shared(A,B,C) private(i,j,k)
    {
        int tt = omp_get_num_threads();
        printf("threads in openmp %d", tt);

#pragma omp for schedule(static)
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++){
                for (j = 0; j < spl; j++) {
                    C[i+j*n] = C[i+j*n] + A[i+k*n] * B[k+j*n];
                }
            }
        }
    }
}

int correct(double *A, double *B, int n) {
    int i, j;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (A[i+j*n] != B[i+j*n]) {
                return 0;
            }
        }
    }
    return 1;
}

void print(double *A, int n) {
    printC(A, n, n);
}

int main (int argc, char* argv[])
{
    int n = 4000;
    n = atoi(argv[1]);
    int master = 0;
    int rank, size;
    int alignment = 16;
    int i;
    // Medidas de tiempo
    // ...............
    double tucpu, tscpu;
    double td1, td2, tr1, tr2; // tiempo despliege y recogida
    double tp1, tp2;	// calculo local, tiempo multiplicacion
    double tot1, tot2;
    // ...............
    // ...............
    MPI_Init (&argc, &argv);      /* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
    printf("Maquina en rango %d de %d\n", rank, size);
    if (rank == master) {
        printf("Matrix size %i \n",n);
        int spl = n/size;
        printf("Sep %i \n",spl);
        // create two matrices of size n
        double *A = (double *) malloc( n*n*sizeof(double) );
        iniMatrix(A,n);
        double *B = (double *) malloc( n*n*sizeof(double) );
        iniMatrix(B,n);
        double *D = (double *) malloc( n*n*sizeof(double) );
        double *C = (double *) malloc( n*n*sizeof(double) );
        // send setting values

        // ctimer(&tot1, &tucpu, &tscpu);

        ctimer(&td1, &tucpu, &tscpu);

        for (i = 0; i < size; i++){
            if (i != master) {
                MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&spl, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&A[0], n*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }
        // .....
        double *AB = (double *) malloc( n*spl*sizeof(double) );
        int indx = 0;
        for (i=0; i<n; i++){
            int dest = i/spl;
            if (dest == master) {
                // Parte que se queda el master
                memcpy(&AB[indx*n], &B[i*n], n*sizeof(double));
                indx++;
            } else {
                // Parte que se reparte a los nodos
                MPI_Send(&B[i*n], n, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            }
        }

        ctimer(&td2, &tucpu, &tscpu);

        double *R  = (double *) malloc( n*spl*sizeof(double) );

        ctimer(&tp1, &tucpu, &tscpu);

        vectMatrix(A, AB, R, n, spl);

        ctimer(&tp2, &tucpu, &tscpu);

        // La parte local a la matriz resultado
        memcpy(&C[0], &R[0], n*spl*sizeof(double));

        ctimer(&tr1, &tucpu, &tscpu);

        // wait for the others to arrive
        for (i = 1; i < size; i++){
            double *S = (double *) malloc( n*spl*sizeof(double) );
            MPI_Recv(S, n*spl, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            memcpy(&C[i*n*spl], &S[0], n*spl*sizeof(double));
        }

        // ctimer(&tot2, &tucpu, &tscpu);

        ctimer(&tr2, &tucpu, &tscpu);

        float tScatter = (float) (td2-td1);
        float tMult    = (float) (tp2-tp1);
        float tGather  = (float) (tr2-tr1);
        float tTotal   = tScatter + tMult + tGather;

        printf("Tiempo scatter %f segundos \n",tScatter);
        printf("Tiempo mult en %d es %f segundos \n",rank, tMult);
        printf("Tiempo gather %f segundos \n",tGather);
        printf("TOTAL %f segundos \n",tTotal);

    } else {
        int n, spl;
        // broadcast...
        MPI_Recv(&n, 1, MPI_INT, master, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&spl, 1, MPI_INT, master, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        double *A = (double *) malloc( n*n*sizeof(double) );
        double *AB = (double *) malloc( n*spl*sizeof(double) );
        double *R  = (double *) malloc( n*spl*sizeof(double) );
        MPI_Recv(A, n*n, MPI_DOUBLE, master, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        int j;
        for (j=0; j<spl; j++){
            MPI_Recv(&AB[j*n], n, MPI_DOUBLE, master, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        ctimer(&tp1, &tucpu, &tscpu);

        vectMatrix(A, AB, R, n, spl);

        ctimer(&tp2, &tucpu, &tscpu);

        MPI_Send(&R[0], n*spl, MPI_DOUBLE, master, 0, MPI_COMM_WORLD);

        printf("Tiempo mult en %d es %f segundos \n",rank, (float) (tp2-tp1));

    }
    // printf( "Hello world from process %d of %d\n", rank, size );
    MPI_Finalize();

    return 0;
}