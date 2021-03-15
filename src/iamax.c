#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>
#include <float.h>


CBLAS_INDEX mnblas_isamax(const int N, const float  *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    float max = FLT_MIN;
    CBLAS_INDEX index = 0;
    for(int i=0; i<N; i+=incX){
        if(float_abs(X[i])>max){
            index = i;
            max=float_abs(X[i]);
        }
    }
    return index;
}

CBLAS_INDEX mnblas_idamax(const int N, const double *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    double max = DBL_MIN;
    CBLAS_INDEX index = 0;
    for(int i=0; i<N; i+=incX){
        if(double_abs(X[i])>max){
            index = i;
            max=double_abs(X[i]);
        }
    }
    return index;
}

CBLAS_INDEX mnblas_icamax(const int N, const void   *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    complexe_float_t* Xc= (complexe_float_t*) X; 
    float max = FLT_MIN;
    CBLAS_INDEX index = 0;
    for(int i=0; i<N; i+=incX){
        float value = float_abs(Xc[i].imaginary+Xc[i].real);
        if(value>max){
            index = i;
            max=value;
        }
    }
    return index;  
}

CBLAS_INDEX mnblas_izamax(const int N, const void   *X, const int incX){
    if(N<0 || incX<0)
        return 0;
    complexe_double_t* Xc= (complexe_double_t*) X; 
    double max = DBL_MIN;
    CBLAS_INDEX index= 0;
    for(int i=0; i<N; i+=incX){
        double value = double_abs(Xc[i].real + Xc[i].imaginary);
        if(value>max){
            index = i;
            max=value;
        }
    }
    return index;  
}