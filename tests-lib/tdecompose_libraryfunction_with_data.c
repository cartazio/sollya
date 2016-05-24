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

int myownlog(mpfi_t result, mpfi_t x, int n, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;
  
  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;
  
  if(n==0) {
    /* Implementation of the neperian logarithm */

    if(mpfi_nan_p(x)) {
      mpfr_t a;
      mpfr_init2(a,mpfi_get_prec(result));
      mpfi_interv_fr(result,a,a);
      mpfr_clear(a);
      return 0;
    }
    if(mpfi_has_zero(x)) {
      mpfr_t a,b;
      mpfr_init2(a,mpfi_get_prec(result));
      mpfr_init2(b,mpfi_get_prec(result));

      mpfr_set_inf(a,-1);
      mpfi_get_right(b,x);
      mpfr_log(b,b,GMP_RNDU);
      mpfi_interv_fr(result, a, b);

      mpfr_clear(a);
      mpfr_clear(b);
      return 0;
    }
    else {
      mpfr_t a,b;
      mpfr_init2(a,mpfi_get_prec(result));
      mpfr_init2(b,mpfi_get_prec(result));

      mpfi_get_left(a,x);
      mpfr_log(a,a, GMP_RNDD);
      mpfi_get_right(b,x);
      mpfr_log(b,b,GMP_RNDU);
      mpfi_interv_fr(result, a, b);

      mpfr_clear(a);
      mpfr_clear(b);
      return 0;
    }
  }


  if(n==1) {
    /* Implementation of 1/x */
    mpfi_inv(result,x);
    return 0;
  }

  if(n==2) {
    /* Implementation of -1/x^2 */
    mpfi_t temp;
    mpfi_init2(temp, mpfi_get_prec(result));

    mpfi_sqr(temp, x);
    mpfi_inv(temp, temp);
    mpfi_neg(result, temp);

    mpfi_clear(temp);
    return 0;
  }

  /* else */
  mpfr_t a,b;
  mpfr_init2(a, mpfi_get_prec(result));
  mpfr_init2(b, mpfi_get_prec(result));

  mpfr_set_inf(a,-1);
  mpfr_set_inf(b, 1);
  mpfi_interv_fr(result,a,b);

  mpfr_clear(a);
  mpfr_clear(b);
  return 0;
}


int myownlog2(mpfi_t result, mpfi_t x, int n) {
  if(n==0) {
    /* Implementation of the neperian logarithm */

    if(mpfi_nan_p(x)) {
      mpfr_t a;
      mpfr_init2(a,mpfi_get_prec(result));
      mpfi_interv_fr(result,a,a);
      mpfr_clear(a);
      return 0;
    }
    if(mpfi_has_zero(x)) {
      mpfr_t a,b;
      mpfr_init2(a,mpfi_get_prec(result));
      mpfr_init2(b,mpfi_get_prec(result));

      mpfr_set_inf(a,-1);
      mpfi_get_right(b,x);
      mpfr_log(b,b,GMP_RNDU);
      mpfi_interv_fr(result, a, b);

      mpfr_clear(a);
      mpfr_clear(b);
      return 0;
    }
    else {
      mpfr_t a,b;
      mpfr_init2(a,mpfi_get_prec(result));
      mpfr_init2(b,mpfi_get_prec(result));

      mpfi_get_left(a,x);
      mpfr_log(a,a, GMP_RNDD);
      mpfi_get_right(b,x);
      mpfr_log(b,b,GMP_RNDU);
      mpfi_interv_fr(result, a, b);

      mpfr_clear(a);
      mpfr_clear(b);
      return 0;
    }
  }


  if(n==1) {
    /* Implementation of 1/x */
    mpfi_inv(result,x);
    return 0;
  }

  if(n==2) {
    /* Implementation of -1/x^2 */
    mpfi_t temp;
    mpfi_init2(temp, mpfi_get_prec(result));

    mpfi_sqr(temp, x);
    mpfi_inv(temp, temp);
    mpfi_neg(result, temp);

    mpfi_clear(temp);
    return 0;
  }

  /* else */
  mpfr_t a,b;
  mpfr_init2(a, mpfi_get_prec(result));
  mpfr_init2(b, mpfi_get_prec(result));

  mpfr_set_inf(a,-1);
  mpfr_set_inf(b, 1);
  mpfi_interv_fr(result,a,b);

  mpfr_clear(a);
  mpfr_clear(b);
  return 0;
}

int main(void) {
  sollya_obj_t g, f, h;
  mpfr_t x,y;
  data_t data = { "Hello world", 0 };
  int (*func)(mpfi_t, mpfi_t, int, void *);
  int deriv;
  void *ptr;
  void (*resDealloc)(void *);

  sollya_lib_init();
  g = sollya_lib_parse_string("sin(_x_^2)^2 * 1/3");
  f = sollya_lib_libraryfunction_with_data(g, "superfunc", myownlog, &data, dealloc);
  sollya_lib_printf("%b (expecting: superfunc)\n", f);
  mpfr_init2(x, 30);
  mpfr_init2(y, 50);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of log(sin(4)^2 * 1/3))\n", y);
  func = NULL; deriv = -1; h = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryfunction_with_data(&func, &deriv, &h, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("h = %b\n", h);
    sollya_lib_printf("derivation order = %d\n", deriv);
    sollya_lib_printf("The function pointer is %s\n", ((func == myownlog) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((ptr == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
    sollya_lib_clear_obj(h);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);
  sollya_lib_clear_obj(g);

  f = sollya_lib_parse_string("diff(superfunc(sin(_x_)))");
  sollya_lib_printf("%b (expecting: (diff(superfunc))(sin(_x_)) * cos(_x_))\n", f);
  sollya_lib_clear_obj(f);
  f = sollya_lib_parse_string("diff(superfunc(_x_))");
  sollya_lib_printf("%b (expecting: diff(superfunc))\n", f);
  mpfr_set_ui(x, 4, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of 1/4)\n", y);
  func = NULL; deriv = -1; h = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryfunction_with_data(&func, &deriv, &h, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("h = %b\n", h);
    sollya_lib_printf("derivation order = %d\n", deriv);
    sollya_lib_printf("The function pointer is %s\n", ((func == myownlog) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((ptr == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
    sollya_lib_clear_obj(h);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);

  f = sollya_lib_parse_string("sin(exp(_x_)) + 3");
  sollya_lib_printf("%b (expecting: sin(exp(_x_)) + 3\n", f);
  func = NULL; deriv = -1; h = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryfunction_with_data(&func, &deriv, &h, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("h = %b\n", h);
    sollya_lib_printf("derivation order = %d\n", deriv);
    sollya_lib_clear_obj(h);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);

  g = sollya_lib_parse_string("sin(_x_^2)^2 * 1/3");
  f = sollya_lib_libraryfunction(g, "superfunc2", myownlog2);
  sollya_lib_printf("%b (expecting: superfunc2(sin(_x_^2)^2 / 3))\n", f);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of log(sin(4)^2 * 1/3))\n", y);
  func = NULL; deriv = -1; h = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryfunction_with_data(&func, &deriv, &h, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("h = %b\n", h);
    sollya_lib_printf("derivation order = %d\n", deriv);
    sollya_lib_clear_obj(h);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);
  sollya_lib_clear_obj(g);

  g = sollya_lib_parse_string("sin(_x_^2)^2 * 1/3");
  f = sollya_lib_libraryfunction_with_data(g, "superfunc3", myownlog, &data, NULL);
  sollya_lib_printf("%b (expecting: superfunc3(sin(_x_^2)^2 / 3))\n", f);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of log(sin(4)^2 * 1/3))\n", y);
  func = NULL; deriv = -1; h = NULL; ptr = NULL; resDealloc = NULL;
  if (sollya_lib_decompose_libraryfunction_with_data(&func, &deriv, &h, &ptr, &resDealloc, f)) {
    sollya_lib_printf("Decomposition of %b succeeded\n", f);
    sollya_lib_printf("h = %b\n", h);
    sollya_lib_printf("derivation order = %d\n", deriv);
    sollya_lib_printf("The function pointer is %s\n", ((func == myownlog) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((ptr == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == NULL) ? "okay" : "wrong"));
    sollya_lib_clear_obj(h);
  } else {
    sollya_lib_printf("Could not decompose %b\n", f);
  }
  sollya_lib_clear_obj(f);
  sollya_lib_clear_obj(g);

  mpfr_clear(x);
  mpfr_clear(y);
  sollya_lib_close();
  return 0;
}
