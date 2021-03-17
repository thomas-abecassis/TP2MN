#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

int main(){

    float fv1[2] = {2,2};
    

    printf("snrm2 pour {2,2} : %f \n", mnblas_snrm2(2, fv1, 1));

    double fd1[2] = {2,2};

    printf("dnrm2 pour {2,2} : %f \n", mnblas_dnrm2( 2, fd1, 1));

    complexe_float_t cf11 = {1,1};
    complexe_float_t cf12 = {1,1};

    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_1[0]=cf11; vect_cf_1[1]=cf12; 

    printf("scnrm2 pour  {(1,1),(1,1)} : %f \n", mnblas_scnrm2(2, vect_cf_1, 1));

    complexe_double_t cd11 = {-1,1};
    complexe_double_t cd12 = {1,1};

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_1[0]=cd11; vect_cd_1[1]=cd12;

    printf("dznrm2 pour {(-1,1),(1,1)} : %f \n", mnblas_dznrm2(2, vect_cd_1, 1));
    return 0;
}