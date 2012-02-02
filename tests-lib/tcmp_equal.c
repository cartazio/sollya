#include <sollya.h>

int callback(int message) {
  if (message == SOLLYA_MSG_TEST_RELIES_ON_FP_RESULT_FAITHFUL_BUT_TOGGLING)
    sollya_lib_printf("The following test relies on floating point evaluation and is not certified\n");
  else
    sollya_lib_printf("Unexpected warning message %d\n", message);
  return 0;
}

int void_callback(int message) {
  return 0;
}

sollya_obj_t sollya_lib_error() {
  sollya_obj_t a, res;
  int (*previous_callback)(int);
  previous_callback = sollya_lib_get_msg_callback();

  sollya_lib_install_msg_callback(void_callback);
  a = SOLLYA_CONST(1);
  res = sollya_lib_and(a, a);
  sollya_lib_clear_obj(a);
  sollya_lib_uninstall_msg_callback();
  sollya_lib_install_msg_callback(previous_callback);
  return res;
}

int main(void) {
  sollya_obj_t a, b, c, d, res;
  mpfr_t tmp1, tmp2;
  double nan = 0./0.;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback);

  /* Equality between numbers */
  a = SOLLYA_CONST(17);
  b = SOLLYA_CONST(42);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(42);
  b = SOLLYA_CONST(17);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(17);
  b = SOLLYA_CONST(17);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* Equality between arithmetical expressions */
  a = SOLLYA_SIN(SOLLYA_DIV(SOLLYA_PI, SOLLYA_CONST(3)));
  b = SOLLYA_DIV(SOLLYA_SQRT(SOLLYA_CONST(3)), SOLLYA_CONST(2));
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_SIN(SOLLYA_DIV(SOLLYA_PI, SOLLYA_CONST(3)));
  b = SOLLYA_DIV(SOLLYA_SQRT(SOLLYA_CONST(3)), SOLLYA_CONST(2));
  b = SOLLYA_MUL(b, SOLLYA_ADD(SOLLYA_CONST(1), sollya_lib_parse_string("1b-200;")));
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* Equality involving errors */
  a = sollya_lib_error();
  b = SOLLYA_CONST(17);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = sollya_lib_error();
  b = sollya_lib_error();
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(nan);
  b = sollya_lib_error();
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* Equality involving NaN */
  a = SOLLYA_CONST(nan);
  b =  SOLLYA_CONST(nan);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(nan);
  b =  SOLLYA_CONST(1);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* Commutativity */
  a = SOLLYA_ADD(SOLLYA_CONST(2), SOLLYA_X_);
  b = SOLLYA_ADD(SOLLYA_X_, SOLLYA_CONST(2));
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* Equality between equal intervals but with bounds at different precisions */
  a = sollya_lib_parse_string("[1;2];");
  mpfr_init2(tmp1, 20); /* not the current precision */
  mpfr_init2(tmp2, 30); /* not the current precision */
  mpfr_set_ui(tmp1, 1, GMP_RNDN);
  mpfr_set_ui(tmp2, 2, GMP_RNDN);
  b = sollya_lib_range_from_bounds(tmp1, tmp2);
  mpfr_clear(tmp1); mpfr_clear(tmp2);
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b (even when bounds are stored at different precisions)\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* More involved expression */
  a = sollya_lib_parse_string("1 + exp(sin(x + log(x^2)));");
  c = SOLLYA_CONST(40);
  sollya_lib_set_prec(c);
  sollya_lib_clear_obj(c);
  b = sollya_lib_parse_string("sin(x + log(x^2));");
  b = SOLLYA_ADD(SOLLYA_CONST(1), SOLLYA_EXP(b));
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b (even when constructed at different precisions)\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  /* Equality between strings */
  /* Structural equality between strings */
  a = sollya_lib_string("Hello!");
  b = sollya_lib_string("Hella!");
  res = sollya_lib_cmp_equal(a, b);
  sollya_lib_printf("%b == %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = sollya_lib_string("Hello!");
  b = sollya_lib_string("ello!");
  c = sollya_lib_string("H");
  d = sollya_lib_concat(c, b);
  res = sollya_lib_cmp_equal(a, d);
  sollya_lib_printf("%b == %b returns %b\n", a, d, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);
  sollya_lib_clear_obj(d);
  sollya_lib_clear_obj(res);

  sollya_lib_uninstall_msg_callback();
  sollya_lib_close();
  return 0;
}
