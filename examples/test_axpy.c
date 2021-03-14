#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

int main(){

    float fv1[2] = {2,2};
    float fv2[2] = {1,1};
    mnblas_saxpy(2, 2, fv1, 1, fv2, 1);

    printf("axpy pour 2, {2,2} et {1,1} : %f %f \n", fv1[0], fv1[1]);


    return 0;
}