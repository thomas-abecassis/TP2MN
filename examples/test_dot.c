#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;

vfloat vec1, vec2 ;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 float res ;
 int i ;

  float v1[2]={1,1};
  float v2[2]={2,2};

  printf("résultat : %f \n", mncblas_sdot(2, v1, 1, v2, 1));

  complexe_float_t* c = malloc(sizeof(complexe_float_t));


  complexe_float_t* c1 = malloc(sizeof(complexe_float_t));
  c1->real = 2; c1->imaginary=2;

  complexe_float_t** vc = malloc(sizeof(complexe_float_t*) * 2);
  vc[0] = c1; vc[1] = c1;


  mncblas_cdotu_sub(2, vc, 1, vc, 1, c);

  printf("résultat complexe : %f + %fi \n", c->real, c->imaginary);

 init_flop () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     res = 0.0 ;
     
     start = _rdtsc () ;
        res = mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }
}
