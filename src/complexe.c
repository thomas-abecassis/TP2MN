#include "complexe.h"

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

  /* 
     a implementer
  */


  r.real = 0.0 ;
  r.imaginary = 0.0 ;
  
  return r ;
}

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
  {
  complexe_double_t r ;

  /* 
     a implementer
  */
  
  r.real = 0.0 ;
  r.imaginary = 0.0 ;
  
  return r ;
}

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2) {

  complexe_float_t c=mult_complexe_float(c1,conjugue_float(c2));
  c=mult_scalaire_float(c,module_carre_float(c2));
  return c;
}

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2){
  complexe_double_t c=mult_complexe_double(c1,conjugue_double(c2));
  c=mult_scalaire_double(c,module_carre_double(c2));
  return c;
}
  


