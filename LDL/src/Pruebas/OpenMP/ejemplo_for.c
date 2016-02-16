//
// Created by myhay on 4/11/15.
//

/* C Example */
#include <omp.h>
#include <stdio.h>

int main(){

    int iam =0, np = 1, i=0;

    #pragma omp parallel private(iam, np,i)
    {
        #if defined (_OPENMP)
            np = omp_get_num_threads();
            iam = omp_get_thread_num();
        #endif

        printf("Hello from thread %d out of %d\n",iam,np);

        #pragma omp for
        for(i=0;i<(np*2);i++)
        {
            printf("Thread %d, contador %d \n",iam,i);
        }
    }
}