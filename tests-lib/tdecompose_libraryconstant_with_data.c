#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_struct_t {
  char text[32];
  int  counter;
} data_t;

void dealloc(void *data) {
  (void) data;
  return;
}

void euler_gamma(mpfr_t res, mp_prec_t prec, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;
  
  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;
  
  mpfr_set_prec(res, prec);
  mpfr_const_euler(res, GMP_RNDN);
  return;
}

void euler_gamma2(mpfr_t res, mp_prec_t prec) {
  mpfr_set_prec(res, prec);
  mpfr_const_euler(res, GMP_RNDN);
  return;
}

int main(void) {
  sollya_obj_t f;
  mpfr_t x,y;
  data_t data = { "Hello world", 0 };
  void (*func)(mpfr_t, mp_prec_t, void *);
  void *ptr;
  void (*resDealloc)(void *);


  sollya_lib_init();

  f = sollya_lib_libraryconstant_with_data("superconst", euler_gamma, &data, dealloc);
  sollya_lib_printf("%b (expecting: superconst)\n", f);
  mpfr_init2(x, 30);
  mpfr_init2(y, 50);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: 0.5772...)\n", y);
  func = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryconstant_with_data(&func, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("The function pointer is %s\n", ((func == euler_gamma) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((ptr == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);

  f = sollya_lib_parse_string("superconst");
  sollya_lib_printf("%b (expecting: superconst)\n", f);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: 0.5772...)\n", y);
  func = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryconstant_with_data(&func, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("The function pointer is %s\n", ((func == euler_gamma) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((ptr == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);

  f = SOLLYA_PI;
  sollya_lib_printf("%b (expecting: pi)\n", f);
  func = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryconstant_with_data(&func, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);

  f = sollya_lib_libraryconstant("superconst2", euler_gamma2);
  sollya_lib_printf("%b (expecting: superconst2)\n", f);
  func = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryconstant_with_data(&func, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);

  f = sollya_lib_libraryconstant_with_data("superconst3", euler_gamma, &data, NULL);
  sollya_lib_printf("%b (expecting: superconst)\n", f);
  func = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryconstant_with_data(&func, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("The function pointer is %s\n", ((func == euler_gamma) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((ptr == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == NULL) ? "okay" : "wrong"));
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);


  mpfr_clear(x);
  mpfr_clear(y);
  sollya_lib_close();
  return 0;
}
