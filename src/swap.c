#include "mnblas.h"
#include "complexe.h"

void mncblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

}

void mncblas_cswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t* save;
  complexe_float_t** Xc = X;
  complexe_float_t** Yc = Y;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Yc [j] ;
      Yc [j] = Xc [i] ;
      Xc [i] = save ;
    }

}

void mncblas_zswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t* save;
  complexe_double_t** Xc = X;
  complexe_double_t** Yc = Y;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Yc [j] ;
      Yc [j] = Xc [i] ;
      Xc [i] = save ;
    }
}

