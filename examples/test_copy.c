#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

int main(){
    float* vect_float_1 = malloc(sizeof(float) * 2);
    vect_float_1[0]=1; vect_float_1[1]=2;
    float* vect_float_2 = malloc(sizeof(float) * 2);
    mncblas_scopy(2,vect_float_1,1,vect_float_2,1);


    printf("vecteur float 1 : {%f, %f} \n", vect_float_1[0], vect_float_1[1]);
    printf("vecteur float 2 : {%f, %f} \n", vect_float_2[0], vect_float_2[1]);

    double* vect_double_1 = malloc(sizeof(double) * 2);
    vect_double_1[0]=1; vect_double_1[1]=2;
    double* vect_double_2 = malloc(sizeof(double) * 2);
    mncblas_dcopy(2,vect_double_1,1,vect_double_2,1);

    printf("vecteur double 1 : {%f, %f} \n", vect_double_1[0], vect_double_1[1]);
    printf("vecteur double 2 : {%f, %f} \n", vect_double_2[0], vect_double_2[1]); 



    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_1[0].real=1.1;
    vect_cf_1[0].imaginary=1.5;

    vect_cf_1[1].real=2.1;
    vect_cf_1[1].imaginary=2.5;

    complexe_float_t* vect_cf_2 = malloc(sizeof(complexe_float_t) * 2);
    mncblas_ccopy(2,vect_cf_1,1,vect_cf_2,1);

    printf("vecteur float 1 : {%f + %fi, %f + %fi} \n", vect_cf_1[0].real, vect_cf_1[0].imaginary, vect_cf_1[1].real, vect_cf_1[1].imaginary);
    printf("vecteur float 2 : {%f + %fi, %f + %fi} \n", vect_cf_2[0].real, vect_cf_2[0].imaginary, vect_cf_2[1].real, vect_cf_2[1].imaginary);  

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_1[0].real=1.1;
    vect_cd_1[0].imaginary=1.5;

    vect_cd_1[1].real=2.1;
    vect_cd_1[1].imaginary=2.5;

    complexe_double_t* vect_cd_2 = malloc(sizeof(complexe_double_t) * 2);
    mncblas_zcopy(2,vect_cd_1,1,vect_cd_2,1);

    printf("vecteur double 1 : {%f + %fi, %f + %fi} \n", vect_cd_1[0].real, vect_cd_1[0].imaginary, vect_cd_1[1].real, vect_cd_1[1].imaginary);
    printf("vecteur double 2 : {%f + %fi, %f + %fi} \n", vect_cd_2[0].real, vect_cd_2[0].imaginary, vect_cd_2[1].real, vect_cd_2[1].imaginary);  

    return 0;

}