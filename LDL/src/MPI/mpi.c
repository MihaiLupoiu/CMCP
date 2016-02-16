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
#include <mpi.h>

#include <string.h>
#include <unistd.h>

double frobeniusNorm(double *A, double *LDLT, size_t size1, size_t size2)
{
    double result = 0.0;
    int i,j;
    for(i = 0; i < size1; ++i)
    {
        for(j = 0; j < size2; ++j)
        {
            double value = A[i*size2+j] - LDLT[i*size2+j];
            result += value * value;
        }
    }
    return sqrt(result);
}

double norm(double* LDL, double* A, size_t size){

    double * LD = (double *) calloc(size*size, sizeof(double));

    int i,j,k = 0;

    for (i = 0; i < size; ++i) {
        for (j = 0; j <= i; j++) {
            if (i == j) {
                LD[i*size+j] = LDL[i*size+j];
            } else {
                LD[i*size+j] += LDL[i*size+j] * LDL[j*size+j];
            }
        }
    }

    double * LT = (double*)calloc(size*size, sizeof(double));
    for (i = 0; i < size; ++i) {
        for (j = i; j < size; j++) {
            if (i == j) {
                LT[i*size+j] = 1;
            } else {
                LT[i*size+j] = LDL[j*size+i];
            }

        }
    }

    double * LDLT = (double*)calloc(size*size, sizeof(double));
    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; j++) {
            for (k = 0; k < size; k++) {
                LDLT[i*size+j] = LDLT[i*size+j] + LD[i*size+k]*LT[k*size+j];
            }
        }
    }

    double resultado = frobeniusNorm(A,LDLT, size,size);

    free(LD);
    free(LT);
    free(LDLT);

    return resultado;
}

int main(int argc, char *argv[]){

    double time_start, time_end;
    double time_sending_data, time_procesing, time_recovering_data;

    MPI_Status status;
    int pid, nprocs;
    int root = 0;

    MPI_Init (&argc, &argv);      /* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);        /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);        /* get number of processes */

    //printf("Maquina %d de %d\n", pid, nprocs);
    if (pid == root){
        /**************************** master ************************************/

        // ================================ //
        double *matrix;
        double *originalMatrix;
        int matrix_size;
        // ================================ //

        if(argc !=2){
            printf("Enter the size of matrix (N x N) where N = ");
            scanf("%d",&matrix_size);
        }
        else{
            matrix_size=atol(argv[1]);
        }

        matrix = (double *) calloc(matrix_size * matrix_size, sizeof(double));
        originalMatrix = (double *) calloc(matrix_size * matrix_size, sizeof(double));

        initializeCust(matrix, matrix_size);
        copyMatrixCustom(matrix, originalMatrix, matrix_size);

        //printMatrixCustom(matrix,matrix_size,matrix_size);
        //printMatrixCustom(originalMatrix,matrix_size,matrix_size);

        time_start=MPI_Wtime();
        // ================================ //
        int s = matrix_size/nprocs;

        int i,j,k;
        for (i = 1; i < nprocs; ++i) {
            MPI_Send(&s, 1, MPI_INT, i, 0,MPI_COMM_WORLD);
            MPI_Send(&matrix_size, 1, MPI_INT, i, 0,MPI_COMM_WORLD);
        }

        // ================ Declarar A_LOCAL ================ //

        int rows, columns;
        if(matrix_size%nprocs != 0){
            rows = s+1;
        }else{
            rows = s;
        }
        columns = matrix_size;

        //printf("%d, %d, %d\n",pid,(int)rows,(int)columns);

        double *A = (double *)calloc(rows*columns,sizeof(double));

        // ================ Declarar A_LOCAL ================ //

        // ================ Reparticion de MATRIX ================ //

        int i_local = 0;
        for (i=0; i<matrix_size; i++){
            int dest = i%nprocs;
            if (dest == root) {
                //printf(" Local => %f, %f\n",matrix[i*columns+0],matrix[i*columns+1]);
                memcpy(&A[i_local*columns], &matrix[i*columns], columns*sizeof(double));
                i_local++;
            } else {
                //printf(" Envío destino %d => %f, %f\n",dest,matrix[i*columns],matrix[i*columns+1]);
                MPI_Send(&(matrix[i*columns]), columns, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            }
        }
        // ================ Reparticion de Matriz ================ //

        time_end=MPI_Wtime();

        time_sending_data = time_end-time_start;

        //printMatrixCustomCPU(A,rows ,columns,pid);

        // ================ Procesado ================ //

        double *D = (double *)calloc(1*columns,sizeof(double));
        double *v = (double *)calloc(1*columns,sizeof(double));

        time_start=MPI_Wtime();

        for (j = 0; j < matrix_size ; ++j) {
            if(pid == j%nprocs){

                int j_local = j/nprocs;

                for (i = 0; i < j ; ++i) {
                    //printf("A= %f\n",A[i+j_local*columns]);
                    v[i] = D[i]*A[i+j_local*columns]; //Dependencia de datos, necesita tener la D (A[i][i]) y L (A[j][i]) anterior.
                    //V es privada.
                    //printf("v[%d] = %f, D[%d]= %f A[%d][%d] = %f\n\n",i,v[i],i,D[i],j_local,i,A[i+j_local*columns]);

                }

                double ts = 0;
                for (k = 0; k < j ; ++k) {
                    ts = ts + A[k+j_local*columns]*v[k]; // //Dependencia de datos, necesita tener la L (A[j][k]) anterior.
                    //printf("ts = %f, V[%d]= %f A[%d][%d] = %f\n\n",ts,k,v[k],j_local,k,A[k+j_local*columns]);

                }

                v[j] = A[j+j_local*columns]-ts;
                //printf("A [%d][%d] \n",j_local,j);
                //printf("v[%d] = %f, A[%d][%d] = %f, ts= %f \n\n",j,v[j],j,j_local,A[j_local+j*columns],ts);

                A[j+j_local*columns] = v[j]; // Se escribe la nueva D (A[j][j])
            }

            int sender = j%nprocs;

            if(matrix_size > j+1){
                //printf("Envio %d %d | %d, %d\n",pid,sender,j+1,matrix_size);
                MPI_Bcast(&(v[0]), columns, MPI_DOUBLE, sender, MPI_COMM_WORLD);
            }//else{
                //printf("BREAK %d | %d, %d\n",pid,j+1,matrix_size);
            //}

            D[j] = v[j]; // Se escribe la nueva D local
            //printf("D= %f\n",D[j]);


            //printMatrixCustomCPU(A,rows,columns,pid);

            int start, end;

            if(pid <= j%nprocs){
                start = (j/nprocs)+1;
            }else{
                start = (j/nprocs);
            }

            if( (matrix_size%nprocs <= pid) && (matrix_size%nprocs != 0) ){
                end = s-1;
            }else{
                end = s;
            }

            //printf("pid: %d, start = %d, end = %d\n",pid,start,rows);

            int i_local; //i_local < end funciona pero es el bueno
            //for (i_local = start; i_local < end; ++i_local) {
            for (i_local = start; i_local < rows; ++i_local) {
                double ts = 0;
                //printf("i_local = %d, j = %d\n",i_local,j);

                for (k = 0; k < j ; ++k) {
                    //printf("i_local = %d, k = %d, pos = %d \n",i_local,k,i_local+k*columns);
                    //printf("A[%d] = %f , v[%d] = %f\n",i_local+k*columns,A[i_local+k*columns],k,v[k]);
                    ts = ts + A[k+i_local*columns]*v[k]; //Dependencia de datos, necesita tener la L (A[i][k]) anterior.

                }
                //printf("===>A[%d][%d] v[%d] = %f , ts = %f\n",i_local,j,j,v[j], ts);

                //printf("=>A[%d][%d] = %f\n",i_local,j,A[j+i_local*columns]);

                A[j+i_local*columns] = (A[j+i_local*columns]-ts)/v[j]; // // Se escribe la nueva L (A[i][j])

                //printf("Result ===>A[%d] = %f \n",i_local+j*columns,A[j+i_local*columns]);
            }

            //printMatrixCustomCPU(A,rows,columns,pid);
            //MPI_Barrier( MPI_COMM_WORLD );

        }

        // ================ Procesado ================ //

        time_end=MPI_Wtime();

        time_procesing = time_end - time_start;

        //printMatrixCustomCPU(A,rows ,columns,pid);

        // ================ Recepcion de Resultado ================ //

        time_start=MPI_Wtime();

        i_local = 0;
        for (i=0; i<matrix_size; i++){

            int origen = i%nprocs;

            if (origen == root) {
                memcpy(&matrix[i*columns], &A[i_local*columns], columns*sizeof(double));
                i_local++;
            } else {
                MPI_Recv(&(matrix[i*columns]), columns, MPI_DOUBLE, origen, 0, MPI_COMM_WORLD,&status);
            }
            //printf("iteracion = %d\n",i);
        }

        time_end=MPI_Wtime();
        time_recovering_data = time_end-time_start;

        //printMatrixCustom(matrix,matrix_size,matrix_size);

        //printf("Frobenius Norm: %20.20f\n", norm(matrix, originalMatrix, matrix_size));

/*
        printf("\n**********************************\n\n");
        printf("Selected :%s\n","MPI");
        printf("Size of Matrix :%d \n",matrix_size);
        printf("Sendig Matrix: %f seconds\n",time_sending_data);
        printf("DECOMPOSE TIME TAKEN : %f seconds\n",time_procesing);
        printf("Recovering result: %f seconds\n",time_recovering_data);
        printf("Total Time Spent: %f seconds\n",time_sending_data+time_procesing+time_recovering_data);
        printf("\n**********************************\n\n");
*/
        printf("%d;",matrix_size);
        printf("%f;",time_sending_data);
        printf("%f;",time_procesing);
        printf("%f;",time_recovering_data);
        printf("%f;",time_sending_data+time_procesing+time_recovering_data);
        printf("%d\n",nprocs);

        // ================ Recepcion de Resultado ================ //

        free(A);
        free(D);
        free(v);
        free(matrix);
        free(originalMatrix);


    }else{ // RECIVE MATRIX
        // **************************************************************** othres ************************************************************************ //

        int j,i,k;
        int s, matrix_size;
        // ================ Declarar A_LOCAL ================ //

        MPI_Recv(&s, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&matrix_size, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);

        int rows, columns;
        if(matrix_size%nprocs != 0){
            rows = s+1;
        }else{
            rows = s;
        }
        columns = matrix_size;

        //printf("%d, %d, %d\n",pid,(int)rows,(int)columns);

        double *A = (double *)calloc(rows*columns,sizeof(double));
        // ================ Declarar A_LOCAL ================ //

        // ================ Reparticion de A_LOCAL ================ //

        int size_transf;
        if( (matrix_size%nprocs <= pid) ){
            size_transf = s-1;
        }else{
            size_transf = s;
        }

        //printf("Size: %d\n",size_transf);

        for (j=0; j<=size_transf; j++){
            //printf("Esperando Recibir pid %d\n",pid);
            MPI_Recv(&A[j*columns], columns, MPI_DOUBLE, root, 0, MPI_COMM_WORLD,&status);
            //printf("Recibido pid = %d: %f %f\n",pid,A[j*columns],A[j*columns+1]);
        }
        // ================ Reparticion de A_LOCAL ================ //

        //printMatrixCustomCPU(A,rows ,columns,pid);

        // ================ Procesado ================ //

        double *D = (double *)calloc(1*columns,sizeof(double));
        double *v = (double *)calloc(1*columns,sizeof(double));

        for (j = 0; j < matrix_size; ++j) {
            //printf("J: %d\n",j);

            if(pid == j%nprocs) {

                int j_local = j / nprocs;

                for (i = 0; i < j; ++i) {
                    v[i] = D[i] * A[i + j_local * columns]; //Dependencia de datos, necesita tener la D (A[i][i]) y L (A[j][i]) anterior.
                    //V es privada.
                    //printf("v[%d] = %f, D[%d]= %f A[%d][%d] = %f\n\n",i,v[i],i,D[i],i,j_local,A[i + j_local * columns]);
                }

                double ts = 0;
                for (k = 0; k < j; ++k) {
                    ts = ts + A[k + j_local * columns] * v[k]; // //Dependencia de datos, necesita tener la L (A[j][k]) anterior.
                    //printf("ts = %f, V[%d]= %f A[%d][%d] = %f\n\n",ts,k,v[k],k,j_local,A[k+j_local*columns]);

                }
                v[j] = A[j + j_local * columns] - ts;

                //printf("A [%d][%d]",j_local,j);
                //printf("v[%d] = %f, A[%d][%d] = %f, ts= %f \n\n",j,v[j],j,j_local,A[j+j_local*columns],ts);

                A[j + j_local * columns] = v[j]; // Se escribe la nueva D (A[j][j])

            }

            int sender = j%nprocs;

            if(matrix_size > j+1){
                //printf("Envio %d %d | %d, %d\n",pid,sender,j+1,matrix_size);
                MPI_Bcast(&(v[0]), columns, MPI_DOUBLE, sender, MPI_COMM_WORLD);
            }//else{
                //printf("BREAK %d | %d, %d\n",pid,j+1,matrix_size);
            //}

            D[j] =v[j]; // Se escribe la nueva D local
            //printf("D= %f, v[%d] = %f\n",D[j],j, v[j]);

            int start, end;

            if(pid <= j%nprocs){
                start = (j/nprocs)+1;
            }else{
                start = (j/nprocs);
            }

           if( (matrix_size%nprocs <= pid) && (matrix_size%nprocs != 0) ){
                    end = s-1;
           }else{
               end = s;
           }

            //printf("=============================================================\n");
            //printf("pid: %d, start = %d, end = %d\n",pid,start,rows);

            int i_local;
            //for (i_local = start; i_local < end; ++i_local) {
            for (i_local = start; i_local < rows; ++i_local) {
                double ts = 0;
                //printf("i_local = %d, j = %d\n",i_local,j);

                for (k = 0; k < j ; ++k) {
                    //printf("i_local = %d, k = %d, pos = %d \n",i_local,k,i_local+k*columns);
                    //printf("A[%d] = %f , v[%d] = %f\n",i_local+k*columns,A[i_local+k*columns],k,v[k]);
                    ts = ts + A[k+i_local*columns]*v[k]; //Dependencia de datos, necesita tener la L (A[i][k]) anterior.

                }
                //printf("===>A[%d][%d] v[%d] = %f , ts = %f\n",i_local,j,j,v[j], ts);
                //printf("=>A[%d][%d] = %f\n",i_local,j,A[j+i_local*columns]);

                A[j+i_local*columns] = (A[j+i_local*columns]-ts)/v[j]; // // Se escribe la nueva L (A[i][j])
                //printf("Result ===>A[%d] = %f \n",i_local+j*columns,A[j+i_local*columns]);
            }

            //printMatrixCustomCPU(A,rows,columns,pid);
            //MPI_Barrier( MPI_COMM_WORLD );
        }

        // ================ Procesado ================ //

        //printMatrixCustomCPU(A,rows ,columns,pid);


        // ================ Envío de Resultado ================ //

        //printf("S: %d\n",s);
        for (j=0; j<=size_transf; j++){
            //printf("Envío: %d, %d\n",pid+(j*s)+j,matrix_size);

            //if( (pid+(j*s)+j) <= matrix_size){
                //printf("%d, %d\n",pid+(j*s)+j,matrix_size);
                //printf("ENVIANDO j=%d \n",j);
                MPI_Send(&A[j*columns], columns, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
            //}
        }
        // ================ Envío de Resultado ================ //

        free(A);
        free(D);
        free(v);
    }

    //MPI_Barrier( MPI_COMM_WORLD );

    MPI_Finalize();

    return 0;
}
