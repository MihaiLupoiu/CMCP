#include <pthread.h>
#include <iostream>
#include <cstdlib>

#define NUM_THREADS 4

typedef struct
{
    float* a;
    float* b;
    unsigned char pr;
    unsigned int n;
}
dotdata_t;

void print_vec(float* v, int n)
{
    std::cout << "[ ";

    for (unsigned int i = 0; i < n; i++)
    {
        std::cout << v[i] << ",";
    }

    std::cout << " ]";
}

void *dot(void *dotData)
{
    dotdata_t *data = (dotdata_t*) dotData;
    float* a = data->a;
    float* b = data->b;
    int pr = data->pr;
    int n = data->n;

    int tb = n/NUM_THREADS;

    for (int i = 0; i < tb; i++)
    {
        int f = i + pr*tb;
        ts[pr] = ts[pr] + a[f] * b[f];
    }

    return NULL;
}

float ts[NUM_THREADS];

int main (int argc, char* argv[])
{
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    dotdata_t data[NUM_THREADS];

    unsigned int n = 10;

    float* b = new float[n];
    float* a = new float[n];

    std::cout << "Generating vector" << std::endl;

    for (unsigned int i = 0; i < n; i++)
    {
        b[i] = 1;
        a[i] = 0;
    }

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    std::cout << "Launching threads..." << std::endl;

    for (unsigned char i = 0; i < NUM_THREADS; i++)
    {
        data[i].pr = i;
        data[i].n = n;
        data[i].a = a;
        data[i].b = b;
        ts[i] = 0;

        pthread_create(&threads[i], &attr, dot, (void*) &data[i]);
    }

    std::cout << "Waiting for join..." << std::endl;
    float s = 0;

    for (unsigned char i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }

    std::cout << "vector a: ";
    print_vec(a, n);

    std::cout << std::endl << "vector b: ";
    print_vec(b, n);

    std::cout << std::endl << "alpha: " << alpha << std::endl;
    std::cout << "a = b + alpha * y" << std::endl;
    std::cout << "vector a: ";

    pthread_exit(NULL);

    delete[] a;
    delete[] b;
}
