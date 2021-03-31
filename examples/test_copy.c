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

    printf("\n\nteste performance\n");

    printf("test float\n");

    free(vect_float_1);
    init_flop();
    vect_float_1=malloc(sizeof(float)*VECSIZE);
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_float_1[i]=1.1;
        }
        free(vect_float_2);
        vect_float_2=malloc(sizeof(complexe_float_t)*VECSIZE);
        start=_rdtsc();
        mncblas_scopy(VECSIZE,vect_float_1,1,vect_float_2,1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }



    printf("test complex float\n");
    free(vect_cf_1);
    init_flop();
    vect_cf_1=malloc(sizeof(complexe_float_t)*VECSIZE);
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cf_1[i].imaginary=1.1;
           vect_cf_1[i].real=2.2;
        }
        free(vect_cf_2);
        vect_cf_2=malloc(sizeof(complexe_float_t)*VECSIZE);
        start=_rdtsc();
        mncblas_ccopy(VECSIZE,vect_cf_1,1,vect_cf_2,1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }
    printf("test complex double\n");
    free(vect_cd_1);
    init_flop();
    vect_cd_1=malloc(sizeof(complexe_double_t)*VECSIZE);
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cd_1[i].imaginary=1.1;
           vect_cd_1[i].real=2.2;
        }
        free(vect_cd_2);
        vect_cd_2=malloc(sizeof(complexe_double_t)*VECSIZE);
        start=_rdtsc();
        mncblas_zcopy(VECSIZE,vect_cd_1,1,vect_cd_2,1);
        end=_rdtsc();
        printf("nombre de cycle:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }



    return 0;

}