#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

/*
float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  double dot = 0.0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t** Xc = (complexe_float_t**) X;
  complexe_float_t** Yc = (complexe_float_t**) Y;
  complexe_float_t* dotuC = dotu;

  for (i = 0 ; i < N ; i += incX)
    {
      dotuC->real = Xc[i]->real * Yc[j]->real;
      dotuC->imaginary = Xc[i]->imaginary * Yc[j]->imaginary;
      j+=incY ;
    }
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t** Xc = (complexe_float_t**) X;
  complexe_float_t** Yc = (complexe_float_t**) Y;
  complexe_float_t* dotcC = dotc;

  for (i = 0 ; i < N ; i += incX)
    {
      dotcC->real = Xc[i]->real * Yc[j]->real;
      dotcC->imaginary = (-Xc[i]->imaginary) * Yc[j]->imaginary;
      j+=incY;
    }
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t** Xc = (complexe_double_t**) X;
  complexe_double_t** Yc = (complexe_double_t**) Y;
  complexe_double_t* dotuC = (complexe_double_t*) dotu;

  for (i = 0 ; i < N ; i += incX)
    {
      dotuC->real = Xc[i]->real * Yc[j]->real;
      dotuC->imaginary = Xc[i]->imaginary * Yc[j]->imaginary;
      j+=incY ;
    }
}
  
void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t** Xc = (complexe_double_t**) X;
  complexe_double_t** Yc = (complexe_double_t**) Y;
  complexe_double_t* dotcC = (complexe_double_t*) dotc;

  for (i = 0 ; i < N ; i += incX)
    {
      dotcC->real = Xc[i]->real * Yc[j]->real;
      dotcC->imaginary = (-Xc[i]->imaginary) * Yc[j]->imaginary;
      j+=incY ;
    }
}




