#include <sollya.h>

int main(void) {
  sollya_obj_t f, g, h, dom, delta;
  mpfr_t x, y;

  sollya_lib_init();

  f = sollya_lib_parse_string("exp(_x_)");
  g = sollya_lib_parse_string("function(proc(X,n,p) { var res, oldPrec;\n"
                              "                       oldPrec = prec;\n"
                              "                       prec = p!;\n"
                              "                       \"Using procedure function exponential\";\n"
			      "                       res = exp(X);\n"
			      "                       prec = oldPrec!;\n"
			      "                       return res;\n"
			      "                     })");
  dom = sollya_lib_parse_string("[-5;5]");
  delta = sollya_lib_parse_string("[0]");

  h = sollya_lib_annotatefunction(f, g, dom, delta, NULL);

  mpfr_init2(x, 24);
  mpfr_init2(y, 24);

  mpfr_set_si(x, 2, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, h, x, NULL);

  mpfr_set_si(x, 6, GMP_RNDN);
  sollya_lib_evaluate_function_at_point(y, h, x, NULL);

  mpfr_clear(y);
  mpfr_clear(x);

  sollya_lib_clear_obj(f);
  sollya_lib_clear_obj(g);
  sollya_lib_clear_obj(dom);
  sollya_lib_clear_obj(delta);
  sollya_lib_clear_obj(h);
  
  sollya_lib_close();
  return 0;
}
