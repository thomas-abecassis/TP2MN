#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>
#include <math.h>


void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                 const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY){

    int row = N;
    int col = M;

    for (int i= 0; i< row*col ; i += incX){
        float value = 0;
        int t = col;
        if(TransA != MNCblasNoTrans)
            t = row;
        for (int k=0; k<t; k++){
            if(TransA==MNCblasNoTrans)
                value+=alpha*A[i*col+k]*X[i];
            else
                value+=alpha*A[i+k*col]*X[i];
        }
        Y[i]=value+beta*Y[i];
    }
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){
    int row = N;
    int col = M;

    for (int i= 0; i< row*col ; i += incX){
        double value = 0;
        int t = col;
        if(TransA != MNCblasNoTrans)
            t = row;
        for (int k=0; k<t; k++){
            if(TransA==MNCblasNoTrans)
                value+=alpha*A[i*col+k]*X[i];
            else
                value+=alpha*A[i+k*col]*X[i];
        }
        Y[i]=value+beta*Y[i];
    }
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY){

    int row = N;
    int col = M;
    complexe_float_t* Xc= (complexe_float_t*) X;
    complexe_float_t* Yc= (complexe_float_t*) Y;
    complexe_float_t* Ac= (complexe_float_t*) A;
    complexe_float_t Alphac= *(complexe_float_t*) alpha;
    complexe_float_t Betac= *(complexe_float_t*) beta;

    for (int i= 0; i< row ; i += incX){
        complexe_float_t value = {0,0};
        int t = col;
        if(TransA != MNCblasNoTrans)
            t = row;
        for (int k=0; k<t; k++){
            if(TransA==MNCblasNoTrans){
                value=add_complexe_float(value,mult_complexe_float(mult_complexe_float(Ac[i*col+k],Alphac), Xc[i]));
            }
            else if(TransA==MNCblasTrans)
                value=add_complexe_float(value, mult_complexe_float(mult_complexe_float(Ac[i+k*col],Alphac), Xc[i]));
            else
                value=add_complexe_float(value,mult_complexe_float(mult_complexe_float(conjugue_float(Ac[i+k*col]),Alphac), Xc[i]));
        }
        Yc[i]=add_complexe_float(value, mult_complexe_float(Yc[i],Betac));
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY){

    int row = N;
    int col = M;
    complexe_double_t* Xc= (complexe_double_t*) X;
    complexe_double_t* Yc= (complexe_double_t*) Y;
    complexe_double_t* Ac= (complexe_double_t*) A;
    complexe_double_t Alphac= *(complexe_double_t*) alpha;
    complexe_double_t Betac= *(complexe_double_t*) beta;
    
    for (int i= 0; i< row*col ; i += incX){
        complexe_double_t value = {0,0};
        int t = col;
        if(TransA != MNCblasNoTrans)
            t = row;
        for (int k=0; k<t; k++){
            if(TransA==MNCblasNoTrans)
                value=add_complexe_double(value, mult_complexe_double(mult_complexe_double(Ac[i*col+k],Alphac), Xc[i]));
            else if(TransA==MNCblasTrans)
                value=add_complexe_double(value, mult_complexe_double(mult_complexe_double(Ac[i+k*col],Alphac), Xc[i]));
            else
                value=add_complexe_double(value, mult_complexe_double(mult_complexe_double(conjugue_double(Ac[i+k*col]),Alphac), Xc[i]));
        }
        Yc[i]=add_complexe_double(value, mult_complexe_double(Yc[i],Betac));
    }
}