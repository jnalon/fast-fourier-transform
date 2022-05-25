/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is the implementation file for the testing and time measuring functions.
 **************************************************************************************************/

// Include necessary libraries:
#include "test_it.h"


void test_it(void (*f)(Complex *, Complex *, int), int size)
{
    Complex x[1024], X[1024];                          // Vectors will at most 1024 samples;

    for(int i=0; i<size; i++) {                        // Initialize vectors;
        x[i] = cmplx(i, 0);
        X[i] = cmplx(0, 0);
    }
    f(x, X, size);                                     // Perform the transform;
    printf("N = %d | Input | Output:\n", size);        // Show the results;
    for(int i=0; i<size; i++)
        printf("  %2d: (%8.4f, %8.4f) | (%8.4f, %8.4f)\n",
               i, x[i].r, x[i].i, X[i].r, X[i].i);
    printf("------------------------------\n");
}


void float_test_it(void (*f)(float complex *, float complex *, int), int size)
{
    float complex x[1024], X[1024];                    // Vectors will at most 1024 samples;

    for(int i=0; i<size; i++) {                        // Initialize vectors;
        x[i] = i;
        X[i] = 0;
    }
    f(x, X, size);                                     // Perform the transform;
    printf("N = %d | Input | Output:\n", size);        // Show the results;
    for(int i=0; i<size; i++)
        printf("  %2d: (%8.4f, %8.4f) | (%8.4f, %8.4f)\n",
               i, creal(x[i]), cimag(x[i]), creal(X[i]), cimag(X[i]));
    printf("------------------------------\n");
}


void double_test_it(void (*f)(double complex *, double complex *, int), int size)
{
    double complex x[1024], X[1024];                   // Vectors will at most 1024 samples;

    for(int i=0; i<size; i++) {                        // Initialize vectors;
        x[i] = i;
        X[i] = 0;
    }
    f(x, X, size);                                     // Perform the transform;
    printf("N = %d | Input | Output:\n", size);        // Show the results;
    for(int i=0; i<size; i++)
        printf("  %2d: (%8.4f, %8.4f) | (%8.4f, %8.4f)\n",
               i, creal(x[i]), cimag(x[i]), creal(X[i]), cimag(X[i]));
    printf("------------------------------\n");
}


float time_it(void (*f)(Complex *, Complex *, int), int size, int repeats)
{
    Complex x[1024], X[1024];                  // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = cmplx(j, 0);
    t0 = clock();                              // Start counting time;
    for(int j=0; j<repeats; j++)
        (*f)(x, X, size);
    t1 = (float) (clock() - t0) / CLOCKS_PER_SEC / (float) repeats;
    return t1;
}


float float_time_it(void (*f)(float complex *, float complex *, int), int size, int repeats)
{
    float complex x[1024], X[1024];            // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = j;
    t0 = clock();                              // Start counting time;
    for(int j=0; j<repeats; j++)
        (*f)(x, X, size);
    t1 = (float) (clock() - t0) / CLOCKS_PER_SEC / (float) repeats;
    return t1;
}


float double_time_it(void (*f)(double complex *, double complex *, int), int size, int repeats)
{
    double complex x[1024], X[1024];           // Vectors will be at most 1024 samples;
    clock_t t0;                                // Starting time;
    float t1;

    for(int j=0; j<size; j++)                  // Initialize the vector;
        x[j] = j;
    t0 = clock();                              // Start counting time;
    for(int j=0; j<repeats; j++)
        (*f)(x, X, size);
    t1 = (float) (clock() - t0) / CLOCKS_PER_SEC / (float) repeats;
    return t1;
}
