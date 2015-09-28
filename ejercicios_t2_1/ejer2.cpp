#include <pthread.h>
#include <iostream>
#include <cstdlib>

#define NUM_THREADS 4

typedef struct
{
    float* a;
    unsigned char pr;
    unsigned int n;
}
sumdata_t;

float ts[NUM_THREADS];

void *sum(void *sumData)
{
    sumdata_t *data = (sumdata_t*) sumData;

    int tb = data->n/NUM_THREADS;

    for (int i = 0; i < tb; i++)
    {
        int f = i + data->pr*tb;
        ts[data->pr] = ts[data->pr] + data->a[f];
    }

    std::cout << "Suma de pr = " << (int)data->pr << ": " << ts[data->pr] << "\n";
}

int main (int argc, char* argv[])
{
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    sumdata_t data[NUM_THREADS];

    unsigned int n = 500;
    float* a = new float[n];
    srand(time(NULL));

    std::cout << "Generating vector" << std::endl;

    for (unsigned int i = 0; i < n; i++)
    {
        a[i] = 1;
    }

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    std::cout << "Launching threads..." << std::endl;

    for (unsigned char i = 0; i < NUM_THREADS; i++)
    {

        data[i].pr = i;
        data[i].n = n;
        data[i].a = a;
        ts[i] = 0;

        pthread_create(&threads[i], &attr, sum, (void*) &data[i]);
    }

    std::cout << "Waiting for join..." << std::endl;
    float s = 0;

    for (unsigned char i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
        s = s + ts[i];
    }

    std::cout << "Total sum: " << s << std::endl;


    pthread_exit(NULL);
}
