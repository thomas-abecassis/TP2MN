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
printf("test float\n");
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
printf("teste double\n");
   init_flop () ;
 double *vect_double1=malloc(VECSIZE*sizeof(double));
 double *vect_double2=malloc(VECSIZE*sizeof(double));
 for (i = 0 ; i < NB_FOIS; i++)
   {
     for(int j=0;j<VECSIZE;j++){
       vect_double1[j]=1.1;
       vect_double2[j]=2.2;
     }
     res = 0.0 ;
     
     start = _rdtsc () ;
        res = mncblas_ddot (VECSIZE, vect_double1, 1, vect_double2, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }

   printf("test complexe float\n");
    init_flop () ;
    
 complexe_float_t* vect_float_complexe1=malloc(VECSIZE*sizeof(complexe_float_t));
 complexe_float_t* vect_float_complexe2=malloc(VECSIZE*sizeof(complexe_float_t));
 complexe_float_t* test_float=malloc(sizeof(complexe_float_t));
 test_float->imaginary=0;
 test_float->real=0;
 complexe_float_t *testa=malloc(sizeof(complexe_float_t));
 for (i = 0 ; i < NB_FOIS; i++)
   {
     for(int j=0;j<VECSIZE;j++){
       vect_float_complexe1[j].real=1.1;
       vect_float_complexe2[j].real=2.2;

       vect_float_complexe1[j].imaginary=1.1;
       vect_float_complexe2[j].imaginary=2.2;
     }
     res = 0.0 ;
     
     start = _rdtsc () ;
     mncblas_cdotu_sub (VECSIZE, vect_float_complexe1, 1, vect_float_complexe2, 1,testa) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }


 printf("test complex double");
 complexe_double_t* vect_double_complexe1=malloc(VECSIZE*sizeof(complexe_double_t));
 complexe_double_t *vect_double_complexe2=malloc(VECSIZE*sizeof(complexe_double_t));
 complexe_double_t* test_double=malloc(sizeof(complexe_double_t));
 test_double->imaginary=0;
 test_double->real=0;
 complexe_double_t *testd=malloc(sizeof(complexe_double_t));
 for (i = 0 ; i < NB_FOIS; i++)
   {
     for(int j=0;j<VECSIZE;j++){
       vect_double_complexe1[j].real=1.1;
       vect_double_complexe2[j].real=2.2;

       vect_double_complexe1[j].imaginary=1.1;
       vect_double_complexe2[j].imaginary=2.2;
     }
     res = 0.0 ;
     
     start = _rdtsc () ;
     mncblas_zdotc_sub (VECSIZE, vect_double_complexe1, 1, vect_double_complexe2, 1,testd) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }


}
