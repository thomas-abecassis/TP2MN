#include "mnblas.h"
#include "complexe.h"

float float_abs(float x){
    if(x <0)
        return -x;
    return x;
}

double double_abs(double x){
    if(x<0)
        return -x;
    return x;
}

float  mnblas_sasum(const int N, const float *X, const int incX){
    register unsigned int i = 0 ;
    float ret=0;

    for (;i < N  ; i += incX)
    {
        ret += float_abs(X[i]);
    }
    return ret;
}

double mnblas_dasum(const int N, const double *X, const int incX){
    register unsigned int i = 0 ;
    double ret=0;

    for (;i < N  ; i += incX)
    {
        ret += double_abs(X[i]);
    }
    return ret;
}

float  mnblas_scasum(const int N, const void *X, const int incX){
    register unsigned int i = 0 ;
    float ret=0;
    complexe_float_t* Xt = (complexe_float_t*) X; 

    for (;i < N  ; i += incX)
    {
        ret += float_abs(Xt[i].imaginary)+float_abs(Xt[i].real);
    }
    return ret;
}

double mnblas_dzasum(const int N, const void *X, const int incX){
    register unsigned int i = 0 ;
    double ret=0;
    complexe_double_t* Xt = (complexe_double_t*) X; 

    for (;i < N  ; i += incX)
    {
        ret += double_abs(Xt[i].imaginary)+double_abs(Xt[i].real);
    }
    return ret;
}
