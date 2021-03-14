#include "mnblas.h"
#include "complexe.h"

void mnblas_saxpy( const int N, const float a, const float *x, const int incX, float *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
        y[j] += a*x[i];
    }
}

void mnblas_daxpy( const int N, const double a, const double *x, const int incX, double *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
        y[j] += a*x[i];
    }
}

void mnblas_caxpy( const int N, const void* a, const void *x, const int incX, void *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t** Xc = (complexe_float_t**) x;
  complexe_float_t** Yc = (complexe_float_t**) y;
  complexe_float_t alpha= *(complexe_float_t*) a;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      complexe_float_t temp = add_complexe_float(mult_complexe_float(*Xc[i], alpha), *Yc[j]);
      Yc[j] = &temp;
    }
}

void mnblas_zaxpy( const int N, const void* a, const void *x, const int incX, void *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t** Xc = (complexe_double_t**) x;
  complexe_double_t** Yc = (complexe_double_t**) y;
  complexe_double_t alpha = *(complexe_double_t*) a;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      complexe_double_t temp = add_complexe_double(mult_complexe_double(*Xc[i], alpha), *Yc[j]);
      Yc[j] = &temp;
    }
}

