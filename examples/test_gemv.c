#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

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
    return 0;
}