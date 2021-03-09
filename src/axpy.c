#include "mnblas.h"
#include "complexe.h"

void cblas_saxpy( const int N, const float a, const float *x, const int incX, float *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
        y[j] += a*x[i];
    }
}

void cblas_daxpy( const int N, const double a, const float *x, const int incX, double *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
        y[j] += a*x[i];
    }
}

void cblas_caxpy( const int N, const float a, const void *x, const int incX, void *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t** Xc = (complexe_float_t**) X;
  complexe_float_t** Yc = (complexe_float_t**) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
        y[j] = &add_complexe_float(scal_complexe_float(x[i], a), y[j]);
    }
}

void cblas_saxpy( const int N, const float a, const void *x, const int incX, void *y, const int incY){
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t** Xc = (complexe_double_t**) X;
  complexe_double_t** Yc = (complexe_double_t**) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
        y[j] = &add_complexe_double(scal_double_float(x[i], a), y[j]);
    }
}

