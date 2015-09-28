#include <pthread.h>
#include <iostream>
#include <cstdlib>

#define NUM_THREADS 4

typedef struct
{
    float* a;
    float* b;
    float* y;
    float alpha;
    unsigned char pr;
    unsigned int n;
}
saxpydata_t;

void print_vec(float* v, int n)
{
    std::cout << "[ ";

    for (unsigned int i = 0; i < n; i++)
    {
        std::cout << v[i] << ",";
    }

    std::cout << " ]";
}

void *saxpy(void *saxpyData)
{
    saxpydata_t *data = (saxpydata_t*) saxpyData;
    float* a = data->a;

    int tb = data->n/NUM_THREADS;

    for (int i = 0; i < tb; i++)
    {
        int f = i + data->pr*tb;
        a[f] = data->b[f] + data->alpha * data->y[f];
    }

    return NULL;
}

int main (int argc, char* argv[])
{
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    saxpydata_t data[NUM_THREADS];

    unsigned int n = 10;

    float* b = new float[n];
    float* a = new float[n];
    float* y = new float[n];
    float alpha = 2;

    std::cout << "Generating vector" << std::endl;

    for (unsigned int i = 0; i < n; i++)
    {
        b[i] = 1;
        a[i] = 0;
        y[i] = 2;
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
        data[i].y = y;
        data[i].alpha = alpha;

        pthread_create(&threads[i], &attr, saxpy, (void*) &data[i]);
    }

    std::cout << "Waiting for join..." << std::endl;
    float s = 0;

    for (unsigned char i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }

    std::cout << "vector b: ";
    print_vec(b, n);

    std::cout << std::endl << "vector y: ";
    print_vec(y, n);

    std::cout << std::endl << "alpha: " << alpha << std::endl;
    std::cout << "a = b + alpha * y" << std::endl;
    std::cout << "vector a: ";
    print_vec(a, n);
    std::cout << std::endl;

    pthread_exit(NULL);

    delete[] a;
    delete[] y;
    delete[] b;
}
