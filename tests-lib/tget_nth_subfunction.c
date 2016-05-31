#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>

void euler_gamma(mpfr_t res, mp_prec_t prec) {
  mpfr_set_prec(res, prec);
  mpfr_const_euler(res, GMP_RNDN);
  return;
}

int stupid1(mpfi_t result, mpfi_t x, int n) {
  (void)x; /* Avoiding "unused parameter" warning */
  (void)n; /* Avoiding "unused parameter" warning */
  mpfi_set_ui(result, 0);
  return 0;
}

int main(void) {
  sollya_obj_t f, tmp, tmp2, tmp3;
  sollya_obj_t g;
  int n;
  int res;

  sollya_lib_init();

  for(n=0;n<=3;n++) {
    /* Constant */
    f = SOLLYA_PI;
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);


    /* Another constant */
    f = SOLLYA_CONST(3);
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* A library constant constant */
    f = sollya_lib_libraryconstant("superconst", euler_gamma);
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* Free variable */
    f = SOLLYA_X_;
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* Elementary function */
    f = SOLLYA_EXP(SOLLYA_ADD(SOLLYA_SIN(SOLLYA_X_), SOLLYA_CONST(17)));
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* Procedure function */
    tmp = sollya_lib_parse_string("proc(X,n,p) {var res, oldPrec; oldPrec = prec; prec = p!; res = exp(X); prec = oldPrec!; return res; };");
    tmp = sollya_lib_build_function_procedurefunction(SOLLYA_X_, tmp);
    tmp2 = sollya_lib_diff(tmp);
    tmp3 = SOLLYA_SIN(SOLLYA_X_);
    f = sollya_lib_apply(tmp2, tmp3, NULL);
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);
    sollya_lib_clear_obj(tmp);
    sollya_lib_clear_obj(tmp2);
    sollya_lib_clear_obj(tmp3);

    /* Library function */
    f = sollya_lib_build_function_libraryfunction(SOLLYA_COS(SOLLYA_X_), "stupid1", stupid1);
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* arithmetic operator */
    f = SOLLYA_ADD(SOLLYA_X_,
                   sollya_lib_build_function_libraryfunction(SOLLYA_X_, "stupid1", stupid1));
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* A Sollya object that is not an expression */
    f = sollya_lib_absolute();
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* Another one */
    f = sollya_lib_parse_string("[1,2]");
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);

    /* What if f is NULL? */
    f = NULL;
    g = NULL; res = -1;
    res = sollya_lib_get_nth_subfunction(&g, f, n);
    if (res) sollya_lib_printf("Subfunction of %b (n=%d): %b\n", f, n, g);
    else sollya_lib_printf("No subfunction of %b such that n=%d\n", f, n);
    sollya_lib_clear_obj(f);
    if (g != NULL)  sollya_lib_clear_obj(g);
  }

  sollya_lib_close();
  return 0;
}
