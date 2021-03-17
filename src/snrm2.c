#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>
#include <math.h>

float  mnblas_snrm2(const int N, const float *X, const int incX){
    float value = 0;
    for(int i=0; i<N; i+=incX){
        value+=X[i]*X[i];
    }
    return sqrtf(value);
}

double mnblas_dnrm2(const int N, const double *X, const int incX){
    double value = 0;
    for(int i=0; i<N; i+=incX){
        value+=X[i]*X[i];
    }
    return sqrt(value);    
}

float  mnblas_scnrm2(const int N, const void *X, const int incX){
    float value = 0;
    complexe_float_t* Xc= (complexe_float_t*) X;
    for(int i=0; i<N; i+=incX){
        value+=sqrtf(module_carre_float(Xc[i]));
    }
    return sqrtf(value);   
}

double mnblas_dznrm2(const int N, const void *X, const int incX){
    double value = 0;
    complexe_double_t* Xc= (complexe_double_t*) X; 
    for(int i=0; i<N; i+=incX){
        value+=sqrt(module_carre_double(Xc[i]));
    }
    return sqrt(value);    
}