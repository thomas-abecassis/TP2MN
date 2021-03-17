#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

int main(){
    unsigned long long start, end ;
    float res ;

    float fv1[2] = {2,2};
    

    printf("sasum pour {2,2} : %f \n", mnblas_sasum(2, fv1, 1));

    double fd1[2] = {2,2};
    mnblas_dasum( 2, fd1, 1);

    printf("dasum pour {2,2} : %f \n", mnblas_dasum( 2, fd1, 1));

    complexe_float_t cf11 = {1,1};
    complexe_float_t cf12 = {1,1};

    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_1[0]=cf11; vect_cf_1[1]=cf12; 

    printf("scasum pour  {(1,1),(1,1)} : %f \n", mnblas_scasum(2, vect_cf_1, 1));

    complexe_double_t cd11 = {-1,1};
    complexe_double_t cd12 = {1,1};

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_1[0]=cd11; vect_cd_1[1]=cd12;

    printf("dzasum pour {(-1,1),(1,1)} : %f \n", mnblas_dzasum(2, vect_cd_1, 1));
    free(vect_cf_1);
    free(vect_cd_1);


    printf("\n\nteste performance\n");
    printf("test complex float\n");
    vect_cf_1=malloc(sizeof(complexe_float_t)*VECSIZE);
    init_flop();
    
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cf_1[i].imaginary=1.1;
           vect_cf_1[i].real=2.2;
        }
        start=_rdtsc();
        mnblas_scasum(2,vect_cf_1,1);
        end=_rdtsc();
        printf("temps d'execution:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }

    printf("\ntest complex double\n");
    init_flop();
    vect_cd_1=malloc(sizeof(complexe_double_t)*VECSIZE);
    for(int i=0;i<NB_FOIS;i++){
        for(int i=0;i<VECSIZE;i++){
           vect_cd_1[i].imaginary=1.1;
           vect_cd_1[i].real=2.2;
        }
        start=_rdtsc();
        mnblas_dcasum(2,vect_cd_1,1);
        end=_rdtsc();
        printf("temps d'execution:%Ld ",end-start);
        calcul_flop("sdot ", 2 * VECSIZE, end-start);
    }
    return 0;
}