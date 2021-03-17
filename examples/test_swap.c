#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

int main(){
    unsigned long long start, end ;

    float* vect_float_1 = malloc(sizeof(float) * 2);
    vect_float_1[0]=1; vect_float_1[1]=1;
    float* vect_float_2 = malloc(sizeof(float) * 2);
    vect_float_2[0]=2; vect_float_2[1]=2;

    mncblas_sswap(2, vect_float_1, 1, vect_float_2, 1);

    printf("vecteur float 1 : {%f, %f} \n", vect_float_1[0], vect_float_1[1]);
    printf("vecteur float 2 : {%f, %f} \n", vect_float_2[0], vect_float_2[1]);

    double* vect_double_1 = malloc(sizeof(double) * 2);
    vect_double_1[0]=1; vect_double_1[1]=1;
    double* vect_double_2 = malloc(sizeof(double) * 2);
    vect_double_2[0]=2; vect_double_2[1]=2;

    mncblas_dswap(2, vect_double_1, 1, vect_double_2, 1);

    printf("vecteur double 1 : {%f, %f} \n", vect_double_1[0], vect_double_1[1]);
    printf("vecteur double 2 : {%f, %f} \n", vect_double_2[0], vect_double_2[1]);  

    complexe_float_t cf11 = {1,1};
    complexe_float_t cf12 = {1.1, 1.1};

    complexe_float_t cf21 = {2,2};
    complexe_float_t cf22 = {2.2,2.2};

    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_1[0]=cf11; vect_cf_1[1]=cf12;
    complexe_float_t* vect_cf_2 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_2[0]=cf21; vect_cf_2[1]=cf22;    

    mncblas_cswap(2,vect_cf_1, 1, vect_cf_2, 1);

    printf("vecteur float 1 : {%f + %fi, %f + %fi} \n", vect_cf_1[0].real, vect_cf_1[0].imaginary, vect_cf_1[1].real, vect_cf_1[1].imaginary);
    printf("vecteur float 2 : {%f + %fi, %f + %fi} \n", vect_cf_2[0].real, vect_cf_2[0].imaginary, vect_cf_2[1].real, vect_cf_2[1].imaginary);  

    complexe_double_t cd11 = {1,1};
    complexe_double_t cd12 = {1.1, 1.1};

    complexe_double_t cd21 = {2,2};
    complexe_double_t cd22 = {2.2,2.2};

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_1[0]=cd11; vect_cd_1[1]=cd12;
    complexe_double_t* vect_cd_2 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_2[0]=cd21; vect_cd_2[1]=cd22;    

    mncblas_cswap(2,vect_cd_1, 1, vect_cd_2, 1);

    printf("vecteur float 1 : {%f + %fi, %f + %fi} \n", vect_cd_1[0].real, vect_cd_1[0].imaginary, vect_cd_1[1].real, vect_cd_1[1].imaginary);
    printf("vecteur double 2 : {%f + %fi, %f + %fi} \n", vect_cd_2[0].real, vect_cd_2[0].imaginary, vect_cd_2[1].real, vect_cd_2[1].imaginary);  

    printf("\n\nteste performance\n");
    printf("test complex float\n");
    free(vect_cf_1);
    free(vect_cf_2);
    init_flop();
    vect_cf_1=malloc(sizeof(complexe_float_t)*VECSIZE);
    vect_cf_2=malloc(sizeof(complexe_float_t)*VECSIZE);
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cf_1[i].imaginary=1.1;
           vect_cf_1[i].real=2.2;
        }
        start=_rdtsc();
        mncblas_cswap(VECSIZE,vect_cf_1,1,vect_cf_2,1);
        end=_rdtsc();
        printf("temps d'execution:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }
    printf("test complex double\n");
    free(vect_cd_1);
    free(vect_cd_2);
    init_flop();
    vect_cf_1=malloc(sizeof(complexe_double_t)*VECSIZE);
    vect_cf_2=malloc(sizeof(complexe_double_t)*VECSIZE);
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cd_1[i].imaginary=1.1;
           vect_cd_1[i].real=2.2;
        }
        start=_rdtsc();
        mncblas_zswap(VECSIZE,vect_cd_1,1,vect_cd_2,1);
        end=_rdtsc();
        printf("temps d'execution:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }

    return 0;
}