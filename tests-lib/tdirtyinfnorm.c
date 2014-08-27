#include <sollya.h>

#define NB_OF_TESTS 7

int callback(sollya_msg_t msg, void *data) {
  (void)data; /* Avoiding "unused parameter" warning */

  int message = sollya_lib_get_msg_id(msg);
  switch(message) {
  case SOLLYA_MSG_EXPRESSION_IS_CONSTANT:
    sollya_lib_printf("Caught message: a certain expression is constant.\n");
    break;
  case SOLLYA_MSG_DOMAIN_IS_NO_CLOSED_INTERVAL_ON_THE_REALS:
    sollya_lib_printf("Caught message: the domain to compute on must be a closed subset of the real numbers.\n");
    break;
  case SOLLYA_MSG_DOMAIN_IS_REDUCED_TO_A_POINT_WILL_SIMPLY_EVAL:
    sollya_lib_printf("Caught message: the domain is reduced to a point.\n");
    break;
  default:
    sollya_lib_printf("Unexpected warning %d.\n", message);
  }
  return 0;
}

int boolean_not_eq(sollya_obj_t a, sollya_obj_t b) {
  sollya_obj_t tmp;
  int res;

  tmp = sollya_lib_cmp_not_equal(a,b);
  res = sollya_lib_is_true(tmp);
  sollya_lib_clear_obj(tmp);
  return res;
}


int main(void) {
  sollya_obj_t a[NB_OF_TESTS], b[NB_OF_TESTS], c[NB_OF_TESTS];
  int i;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback, NULL);

  a[0] = SOLLYA_CONST(17.0);
  b[0] = sollya_lib_parse_string("[-1;1]");

  a[1] = SOLLYA_X_;
  b[1] = sollya_lib_parse_string("[-1;1.5]");

  a[2] = sollya_lib_parse_string("x^2 + x -1");
  b[2] = sollya_lib_parse_string("[-2;1]");

  a[3] = sollya_lib_parse_string("x^2 + x -1");
  b[3] = sollya_lib_parse_string("[-2;2]");

  a[4] = sollya_lib_parse_string("exp(x)");
  b[4] = sollya_lib_parse_string("[-infty;2]");

  a[5] = sollya_lib_parse_string("erf(x)");
  b[5] = sollya_lib_parse_string("[-infty;infty]");

  a[6] = sollya_lib_parse_string("exp(x)");
  b[6] = sollya_lib_parse_string("[2;2]");

  for (i=0;i<NB_OF_TESTS;i++) {
    c[i] = sollya_lib_dirtyinfnorm(a[i],b[i]);
    sollya_lib_printf("The supremum norm of %b over %b gets approximated with dirtyinfnorm by %b\n",a[i],b[i],c[i]);
  }

  for (i=0;i<NB_OF_TESTS;i++) {
    sollya_lib_clear_obj(a[i]);
    sollya_lib_clear_obj(b[i]);
    sollya_lib_clear_obj(c[i]);
  }

  sollya_lib_close();
  return 0;
}

