#include "complexe.h"

float float_abs(float x){
    if(x <0)
        return -x;
    return x;
}

double double_abs(double x){
    if(x<0)
        return -x;
    return x;
}


complexe_double_t conjugue_double(complexe_double_t c){
  complexe_double_t new_c;
  new_c.real=c.real;
  new_c.imaginary=-c.imaginary;
  return new_c;
}

complexe_float_t conjugue_float(complexe_float_t c){
  complexe_float_t new_c;
  new_c.real=c.real;
  new_c.imaginary=-c.imaginary;
  return new_c;
}

double module_carre_double(complexe_double_t c){
  return (c.real*c.real+c.imaginary*c.imaginary);
}
complexe_double_t mult_scalaire_double(complexe_double_t c, double a){
  complexe_double_t new_c;
  new_c.imaginary=c.imaginary*a;
  new_c.real=c.real*a;
  return new_c;
}

float module_carre_float(complexe_float_t c){
  return (c.real*c.real+c.imaginary*c.imaginary);
}
complexe_float_t mult_scalaire_float(complexe_float_t c, float a){
  complexe_float_t new_c;
  new_c.imaginary=c.imaginary*a;
  new_c.real=c.real*a;
  return new_c;
}

complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;
  r.real = c1.real * c2.real - c1.imaginary*c2.imaginary;
  r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;
  
  return r ;
}

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
  {
  complexe_double_t r ;

  r.real = c1.real * c2.real - c1.imaginary*c2.imaginary;
  r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;
  
  return r ;
}

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2) {

  complexe_float_t c=mult_complexe_float(c1,conjugue_float(c2));
  c=mult_scalaire_float(c,1/module_carre_float(c2));
  return c;
}

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2){
  complexe_double_t c=mult_complexe_double(c1,conjugue_double(c2));
  c=mult_scalaire_double(c,1/module_carre_double(c2));
  return c;
}

complexe_float_t scal_complexe_float(const complexe_float_t c, int x){
  complexe_float_t ret;

  ret.real = c.real * x;
  ret.imaginary = c.imaginary * x;

  return ret;
}
  
complexe_double_t scal_complexe_double(const complexe_double_t c, int x){
  complexe_double_t ret;

  ret.real = c.real * x;
  ret.imaginary = c.imaginary * x;

  return ret;
}

