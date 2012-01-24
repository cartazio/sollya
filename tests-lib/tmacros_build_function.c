#include <sollya.h>

int main(void) {
  sollya_obj_t a,b;
  sollya_lib_init();

  a = SOLLYA_SQRT(SOLLYA_X_); b = sollya_lib_build_function_sqrt(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_EXP(SOLLYA_X_); b = sollya_lib_build_function_exp(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_LOG(SOLLYA_X_); b = sollya_lib_build_function_log(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_LOG2(SOLLYA_X_); b = sollya_lib_build_function_log2(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_LOG10(SOLLYA_X_); b = sollya_lib_build_function_log10(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_SIN(SOLLYA_X_); b = sollya_lib_build_function_sin(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_COS(SOLLYA_X_); b = sollya_lib_build_function_cos(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_TAN(SOLLYA_X_); b = sollya_lib_build_function_tan(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ASIN(SOLLYA_X_); b = sollya_lib_build_function_asin(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ACOS(SOLLYA_X_); b = sollya_lib_build_function_acos(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ATAN(SOLLYA_X_); b = sollya_lib_build_function_atan(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_NEG(SOLLYA_X_); b = sollya_lib_build_function_neg(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ABS(SOLLYA_X_); b = sollya_lib_build_function_abs(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_D(SOLLYA_X_); b = sollya_lib_build_function_double(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_SG(SOLLYA_X_); b = sollya_lib_build_function_single(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_QD(SOLLYA_X_); b = sollya_lib_build_function_quad(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

  a = SOLLYA_HP(SOLLYA_X_); b = sollya_lib_build_function_halfprecision(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_DD(SOLLYA_X_); b = sollya_lib_build_function_double_double(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_TD(SOLLYA_X_); b = sollya_lib_build_function_triple_double(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ERF(SOLLYA_X_); b = sollya_lib_build_function_erf(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ERFC(SOLLYA_X_); b = sollya_lib_build_function_erfc(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_LOG1P(SOLLYA_X_); b = sollya_lib_build_function_log1p(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_EXPM1(SOLLYA_X_); b = sollya_lib_build_function_expm1(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_DE(SOLLYA_X_); b = sollya_lib_build_function_doubleextended(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_CEIL(SOLLYA_X_); b = sollya_lib_build_function_ceil(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_FLOOR(SOLLYA_X_); b = sollya_lib_build_function_floor(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_NEARESTINT(SOLLYA_X_); b = sollya_lib_build_function_nearestint(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_SINH(SOLLYA_X_); b = sollya_lib_build_function_sinh(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_COSH(SOLLYA_X_); b = sollya_lib_build_function_cosh(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_TANH(SOLLYA_X_); b = sollya_lib_build_function_tanh(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ASINH(SOLLYA_X_); b = sollya_lib_build_function_asinh(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ACOSH(SOLLYA_X_); b = sollya_lib_build_function_acosh(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_ATANH(SOLLYA_X_); b = sollya_lib_build_function_atanh(SOLLYA_X_);
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

  a = SOLLYA_ADD(SOLLYA_X_, sollya_lib_constant_from_int(1));
  b = sollya_lib_build_function_add((SOLLYA_X_), sollya_lib_constant_from_int(1));
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_SUB(SOLLYA_X_, sollya_lib_constant_from_int(1)); b = sollya_lib_build_function_sub((SOLLYA_X_), ( sollya_lib_constant_from_int(1)));
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_MUL(SOLLYA_X_, sollya_lib_constant_from_int(2)); b = sollya_lib_build_function_mul((SOLLYA_X_), ( sollya_lib_constant_from_int(2)));
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_DIV(SOLLYA_X_, sollya_lib_constant_from_int(2)); b = sollya_lib_build_function_div((SOLLYA_X_), ( sollya_lib_constant_from_int(2)));
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

a = SOLLYA_POW(SOLLYA_X_, sollya_lib_constant_from_int(2)); b = sollya_lib_build_function_pow((SOLLYA_X_), ( sollya_lib_constant_from_int(2)));
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

  a = SOLLYA_PI(); b = sollya_lib_build_function_pi();
  sollya_lib_printf("%b is %b\n",a,b);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);

  sollya_lib_close();
  return 0;
}
