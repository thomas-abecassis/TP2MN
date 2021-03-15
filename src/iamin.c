#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>
#include <float.h>


CBLAS_INDEX mnblas_isamin(const int N, const float  *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    float min = FLT_MAX;
    CBLAS_INDEX index = 0;
    for(int i=0; i<N; i+=incX){
        if(float_abs(X[i])<min){
            index = i;
            min=float_abs(X[i]);
        }
    }
    return index;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    double min = DBL_MAX;
    CBLAS_INDEX index = 0;
    for(int i=0; i<N; i+=incX){
        if(double_abs(X[i])<min){
            index = i;
            min=double_abs(X[i]);
        }
    }
    return index;
}

CBLAS_INDEX mnblas_icamin(const int N, const void   *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    complexe_float_t* Xc= (complexe_float_t*) X; 
    float min = FLT_MAX;
    CBLAS_INDEX index = 0;
    for(int i=0; i<N; i+=incX){
        float value = float_abs(Xc[i].imaginary+Xc[i].real);
        if(value<min){
            index = i;
            min=value;
        }
    }
    return index;  
}

CBLAS_INDEX mnblas_izamin(const int N, const void   *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    complexe_double_t* Xc= (complexe_double_t*) X; 
    double min = DBL_MAX;
    CBLAS_INDEX index= 0;
    for(int i=0; i<N; i+=incX){
        double value = double_abs(Xc[i].real + Xc[i].imaginary);
        if(value<min){
            index = i;
            min=value;
        }
    }
    return index;  
}