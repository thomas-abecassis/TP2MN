#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    256*256

#define NB_FOIS    10

int main(){

    float fv1[4] = {2,3,2,2};
    float fv2[2] = {1,1};
    float fv3[2] = {1,1};
    mncblas_sgemv(101,111, 2, 2, 1, fv1, 1, fv2, 1, 1, fv3, 1);

    printf("sgemv pour 2, {2,2,3,2}, {1,1} et {1,1} : {%f,%f} \n", fv3[0], fv3[1]);

    double fd1[4] = {2,2,2,2};
    double fd2[2] = {1,1};
    double fd3[2] = {1,1};
    mncblas_dgemv(101,111, 2, 2, 1, fd1, 1, fd2, 1, 1, fd3, 1);

    printf("dgemv pour 2, {2,2,2,2}, {1,1} et {1,1} : {%f,%f} \n", fd3[0], fd3[1]);

    complexe_float_t cf11 = {1,1};
    complexe_float_t cf12 = {1,1};

    complexe_float_t cf21 = {1,1};
    complexe_float_t cf22 = {1,1};

    complexe_float_t cm1 = {2,2};    

    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_1[0]=cf11; vect_cf_1[1]=cf12;
    complexe_float_t* vect_cf_2 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_2[0]=cf21; vect_cf_2[1]=cf22; 
    complexe_float_t* mat = malloc(sizeof(complexe_float_t) * 4);
    mat[0]=cm1; mat[1]=cm1; mat[2]=cm1; mat[3]=cm1;  

    mncblas_cgemv(101,111, 2, 2, &cf11, mat, 1, vect_cf_1, 1, &cf11, vect_cf_2, 1);

    printf("cgemv pour {(2,2),(2,2),(2,2),(2,2)}, {(1,1),(1,1)} et {(1,1),(1,1)} : {(%f, %f), (%f, %f)} \n", vect_cf_2[0].real, vect_cf_2[0].imaginary, vect_cf_2[1].real, vect_cf_2[1].imaginary);

    complexe_double_t cd11 = {1,1};
    complexe_double_t cd12 = {1,1};

    complexe_double_t cd21 = {1,1};
    complexe_double_t cd22 = {1,1};

    complexe_double_t cmd1 = {2,2};    

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_1[0]=cd11; vect_cd_1[1]=cd12;
    complexe_double_t* vect_cd_2 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_2[0]=cd21; vect_cd_2[1]=cd22; 
    complexe_double_t* matd = malloc(sizeof(complexe_double_t) * 4);
    matd[0]=cmd1; matd[1]=cmd1; matd[2]=cmd1; matd[3]=cmd1;  

    mncblas_zgemv(101,111, 2, 2, &cd11, matd, 1, vect_cd_1, 1, &cd11, vect_cd_2, 1);

    printf("zgemv pour {(2,2),(2,2),(2,2),(2,2)}, {(1,1),(1,1)} et {(1,1),(1,1)} : {(%f, %f), (%f, %f)} \n", vect_cf_2[0].real, vect_cf_2[0].imaginary, vect_cf_2[1].real, vect_cf_2[1].imaginary);
    
    printf("\n\ndebut test performance\n");
    
    printf("test float");
    unsigned long long start,end;
    float * vect_float=malloc(sizeof(float)*VECSIZE);
    float* vect_float2=malloc(sizeof(float)*VECSIZE);
    float *result_float=malloc(sizeof(float)*VECSIZE);
    for(int i=0;i<VECSIZE;i++){
        result_float[i]=1;
    }
    init_flop();
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_float[i]=1.1;
           vect_float2[i]=2.2;
           
        }
        start=_rdtsc();
        mncblas_sgemv(101,111, 256, 256, 1, vect_float, 1, vect_float2, 1, 1, result_float, 1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }

    printf("test double\n");
    double * vect_double=malloc(sizeof(double)*VECSIZE);
    double* vect_double2=malloc(sizeof(double)*VECSIZE);
    double* result_double=malloc(sizeof(double)*VECSIZE);
    for(int i=0;i<VECSIZE;i++){
        result_double[i]=1;
    }
    init_flop();
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_double[i]=1.1;
           vect_double2[i]=2.2;
        }
        start=_rdtsc();
        mncblas_dgemv(101,111, 256, 256, 1, vect_double, 1, vect_double2, 1, 1, result_double, 1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }

    printf("test complex float\n");
    complexe_float_t scal = {1,0};
    free(vect_cf_1);
    free(vect_cf_2);
    init_flop();
    vect_cf_1=malloc(sizeof(complexe_float_t)*VECSIZE);
    vect_cf_2=malloc(sizeof(complexe_float_t)*VECSIZE);
    complexe_float_t* result_complex_float=malloc(sizeof(complexe_float_t)*VECSIZE);
    for(int i=0;i<VECSIZE;i++){
        result_complex_float[i].imaginary=1.1;
        result_complex_float[i].real=0.5;
    }
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cf_1[i].imaginary=1.1;
           vect_cf_1[i].real=2.2;
           vect_cf_2[i].imaginary=1.1;
           vect_cf_2[i].real=1.1;
        }
        start=_rdtsc();
        mncblas_cgemv(101,111, 256, 256, &scal, vect_cf_1, 1, vect_cf_2, 1, &scal, result_complex_float, 1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }
    printf("test complex double\n");
    complexe_double_t scald = {1,0};
    free(vect_cd_1);
    free(vect_cd_2);
    init_flop();
    vect_cd_1=malloc(sizeof(complexe_double_t)*VECSIZE);
    vect_cd_2=malloc(sizeof(complexe_double_t)*VECSIZE);
    complexe_double_t* result_complex_double=malloc(sizeof(complexe_double_t)*VECSIZE);
    for(int i=0;i<VECSIZE;i++){
        result_complex_double[i].imaginary=1.1;
        result_complex_double[i].real=0.5;
    }

    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cd_1[i].imaginary=1.1;
           vect_cd_1[i].real=2.2;
           vect_cd_2[i].real=1;
           vect_cd_2[i].imaginary=0.5;
        }
        start=_rdtsc();
        mncblas_zgemv(101,111, 256, 256, &scal, vect_cf_1, 1, vect_cf_2, 1, &scald, result_complex_float, 1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }

    return 0;
}