
typedef struct {
  float real ;
  float imaginary ;
} complexe_float_t ;

typedef struct {
  double real ;
  double imaginary ;
} complexe_double_t ;

complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2) ;

complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2) ;

complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2) ;

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2) ;

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2) ;

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2) ;

complexe_float_t scal_complexe_float(const complexe_float_t c, int x);

complexe_double_t scal_complexe_double(const complexe_double_t c, int x);

float module_carre_float(complexe_float_t c);

double module_carre_double(complexe_double_t c);

float float_abs(float x);

double double_abs(double x);
