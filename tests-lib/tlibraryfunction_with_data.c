#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_struct_t {
  char text[32];
  int  counter;
} data_t;

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

int myownlog2(mpfi_t result, mpfi_t x, int n, void *ptr) {
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

void dealloc(void *ptr) {
  data_t *data;
  data = (data_t *) ptr;
  sollya_lib_printf("Deallocation function called for the data pointer (%s<<<>>>%d)\n", data->text, data->counter);
  return;
}

int main(void) {
  sollya_obj_t X, f, g;
  mpfr_t x,y;
  data_t data = { "Hello world", 0 };

  sollya_lib_init();
  X = sollya_lib_free_variable();

  f = sollya_lib_libraryfunction_with_data(X, "superfunc", myownlog, &data, NULL);
  sollya_lib_printf("%b (expecting: superfunc)\n", f);
  mpfr_init2(x, 30);
  mpfr_init2(y, 50);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of log(2))\n", y);
  sollya_lib_clear_obj(f);
  sollya_lib_clear_obj(X);
  f = sollya_lib_parse_string("diff(superfunc(sin(_x_)))");
  sollya_lib_printf("%b (expecting: (diff(superfunc))(sin(_x_)) * cos(_x_))\n", f);
  sollya_lib_clear_obj(f);
  f = sollya_lib_parse_string("diff(superfunc(_x_))");
  sollya_lib_printf("%b (expecting: diff(superfunc))\n", f);
  mpfr_set_ui(x, 4, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of 1/4)\n", y);
  sollya_lib_clear_obj(f);

  X = sollya_lib_free_variable();
  f = sollya_lib_libraryfunction_with_data(X, "superfunc2", myownlog2, &data, dealloc);
  g = SOLLYA_ADD(SOLLYA_CONST(3), sollya_lib_copy_obj(f));

  sollya_lib_printf("%b (expecting: superfunc2)\n", f);
  mpfr_set_ui(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, f, x, NULL);
  sollya_lib_printf("%v (expecting: approximate value of log(2))\n", y);
  sollya_lib_clear_obj(f);
  sollya_lib_printf("Deallocation function must not yet have been called\n");
  sollya_lib_clear_obj(X);
  sollya_lib_clear_obj(g);

  mpfr_clear(x);
  mpfr_clear(y);
  sollya_lib_close();
  return 0;
}
