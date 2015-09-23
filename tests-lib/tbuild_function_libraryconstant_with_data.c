#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_struct_t {
  char text[32];
  int  counter;
} data_t;

void euler_gamma(mpfr_t res, mp_prec_t prec, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;
  
  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;
  
  mpfr_set_prec(res, prec);
  mpfr_const_euler(res, GMP_RNDN);
  return;
}

int main(void) {
  sollya_obj_t f;
  mpfr_t x,y;
  int i;
  data_t data = { "Hello world", 0 };

  sollya_lib_init();

  f = sollya_lib_build_function_libraryconstant_with_data("superconst", euler_gamma, &data);
  sollya_lib_printf("%b (expecting: superconst)\n", f);
  mpfr_init2(x, 30);
  mpfr_init2(y, 50);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: 0.5772...)\n", y);  
  sollya_lib_clear_obj(f);
  f = sollya_lib_parse_string("superconst");
  sollya_lib_printf("%b (expecting: superconst)\n", f);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: 0.5772...)\n", y);  
  sollya_lib_clear_obj(f);
  mpfr_clear(x);
  mpfr_clear(y);
  sollya_lib_close();
  return 0;
}
