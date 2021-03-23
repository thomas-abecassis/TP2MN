#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

int main(){

    float fv1[2] = {1,1};
    float fv2[2] = {1,1};
    float matrice[4] ={2,2,2,2} 

    mncblas_sgemv(101, 111, 2, const int 2,1 , matrice, 2, fv1, 1, 1, fv2, 1);

    printf("sgemv pour {2,2} : %ld \n", mnblas_isamax(5, fv1, 1));

    double fd1[5] = {2,2};

    printf("idamax pour {2,2} : %ld \n", mnblas_idamax( 5, fd1, 1));

    complexe_float_t cf11 = {2,2};
    complexe_float_t cf12 = {1,1};

    complexe_float_t* vect_cf_1 = malloc(sizeof(complexe_float_t) * 2);
    vect_cf_1[0]=cf11; vect_cf_1[1]=cf12; 

    printf("icamax pour  {(2,2),(1,1)} : %ld \n", mnblas_icamax(2, vect_cf_1, 1));

    complexe_double_t cd11 = {2,2};
    complexe_double_t cd12 = {1,1};

    complexe_double_t* vect_cd_1 = malloc(sizeof(complexe_double_t) * 2);
    vect_cd_1[0]=cd11; vect_cd_1[1]=cd12;

    printf("izamax pour {(2,2),(1,1)} : %ld \n", mnblas_izamax(2, vect_cd_1, 1));
    return 0;
}