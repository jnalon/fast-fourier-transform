/**************************************************************************************************
 * Fast Fourier Transform -- C Version
 * This is the implementation file for the time measuring function.
 **************************************************************************************************/

// Include necessary libraries:
#include "time_it.h"


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


float native_complex_time_it(void (*f)(float complex *, float complex *, int), int size, int repeats)
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
