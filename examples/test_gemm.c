#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

int main(){

    float fv1[6] = {1,2,1,2,1,2};
    float fv2[6] = {1,2,1,2,1,2};
    float fv3[4] = {1,1,1,1};
    mncblas_sgemm(101,111,111, 2, 2, 3, 1, fv1, 1, fv2, 1, 1, fv3, 1);

    printf("sgemm pour {1,2,1,2,1,2}, {1,2,1,2,1,2} et {1,1,1,1} : {%f,%f,%f,%f} \n", fv3[0], fv3[1], fv3[2], fv3[3]);

    double fd1[6] = {2,2,2,2,2,2};
    double fd2[6] = {2,2,2,2,2,2};
    double fd3[4] = {0,0,0,0};
    mncblas_dgemm(101,111,111, 2, 2, 3, 1, fd1, 1, fd2, 1, 1, fd3, 1);

    printf("dgemm pour {2,2,2,2,2,2}, {2,2,2,2,2,2} et {0,0,0,0} : {%f,%f,%f,%f} \n", fd3[0], fd3[1], fd3[2], fd3[3]);

    complexe_float_t scal = {1,0};
    complexe_float_t cf11 = {1,0}; 

    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 6);
    complexe_float_t* vect_cf_2 = malloc(sizeof(complexe_float_t) * 6);
    complexe_float_t* mat = malloc(sizeof(complexe_float_t) * 4);
    for(int i=0; i<6; i++){
        vect_cf_1[i]=cf11;
        vect_cf_2[i]=cf11;
    }
    for(int i=0; i<4; i++){
        mat[i]=cf11;
    }

    mncblas_cgemm(101,111,111, 2, 2, 3, &scal, vect_cf_1, 1, vect_cf_2, 1, &scal, mat, 1);

    printf("dgemm : {(%f,%f),(%f,%f),(%f,%f),(%f,%f)} \n", mat[0].real, mat[0].imaginary, mat[1].real, mat[1].imaginary, mat[2].real, mat[2].imaginary, mat[3].real, mat[3].imaginary);
    
    complexe_double_t scald = {1,0};
    complexe_double_t cd11 = {1,0}; 

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 6);
    complexe_double_t* vect_cd_2 = malloc(sizeof(complexe_double_t) * 6);
    complexe_double_t* matd = malloc(sizeof(complexe_double_t) * 4);

    for(int i=0; i<6; i++){
        vect_cd_1[i]=cd11;
        vect_cd_2[i]=cd11;
    }
    for(int i=0; i<4; i++){
        matd[i]=cd11;
    }

    mncblas_zgemm(101,111,111, 2, 2, 3, &scald, vect_cd_1, 1, vect_cd_2, 1, &scald, matd, 1);

    printf("zgemm pour : {(%f,%f),(%f,%f),(%f,%f),(%f,%f)} \n", matd[0].real, matd[0].imaginary, matd[1].real, matd[1].imaginary, matd[2].real, matd[2].imaginary, matd[3].real, matd[3].imaginary);
    
    return 0;

}