#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>
#include <math.h>


void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc){

    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){ 
            float value =0;
            for(int l=0; l<K; l++){
                value+= alpha * A[i*K+l] * B[l*N+j];
            } 
            C[i*N+j]=value + beta * C[i*N+j];
        }
    }
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){

    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){ 
            float value =0;
            for(int l=0; l<K; l++){
                value+= alpha * A[i*K+l] * B[l*N+j];
            } 
            C[i*N+j]=value + beta * C[i*N+j];
        }
    }
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){

    complexe_float_t* Ac= (complexe_float_t*) A;
    complexe_float_t* Bc= (complexe_float_t*) B;
    complexe_float_t* Cc= (complexe_float_t*) C;
    complexe_float_t Alphac= *(complexe_float_t*) alpha;
    complexe_float_t Betac= *(complexe_float_t*) beta;

    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){ 
            complexe_float_t value ={0,0};
            for(int l=0; l<K; l++){
                value = add_complexe_float(value,mult_complexe_float(mult_complexe_float(Ac[i*K+l],Alphac), Bc[l*N+j]));
            } 
            Cc[i*N+j] = add_complexe_float(value, mult_complexe_float(Cc[i*N+j],Betac));
        }
    }
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){

    complexe_double_t* Ac= (complexe_double_t*) A;
    complexe_double_t* Bc= (complexe_double_t*) B;
    complexe_double_t* Cc= (complexe_double_t*) C;
    complexe_double_t Alphac= *(complexe_double_t*) alpha;
    complexe_double_t Betac= *(complexe_double_t*) beta;

    for(int i=0; i<M; i++){
        for(int j=0; j<N; j++){ 
            complexe_double_t value ={0,0};
            for(int l=0; l<K; l++){
                value = add_complexe_double(value,mult_complexe_double(mult_complexe_double(Ac[i*K+l],Alphac), Bc[l*N+j]));
            } 
            Cc[i*N+j] = add_complexe_double(value, mult_complexe_double(Cc[i*N+j],Betac));
        }
    }
}
