#include <sollya.h>
int callback(int message) {
  if (message = SOLLYA_MSG_EXPR_NOT_CORRECTLY_TYPED)
    sollya_lib_printf("The following expression is not correctly typed.\n");
  else
    sollya_lib_printf("Unexpected warning.\n");
  return 0;
}

int main(void) {
  sollya_obj_t a, b, res;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback);

  a = SOLLYA_CONST(1);
  b = sollya_lib_parse_string("[1,2];");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(4);
  b = sollya_lib_parse_string("[1,2];");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = sollya_lib_parse_string("[1,2];");
  b = sollya_lib_parse_string("[3,4];");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = sollya_lib_parse_string("[2,3];");
  b = sollya_lib_parse_string("[1,4];");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = sollya_lib_parse_string("[1,3];");
  b = sollya_lib_parse_string("[2,4];");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = sollya_lib_string("e");
  b = sollya_lib_string("Hello");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(1);
  b = sollya_lib_string("H1llo");
  res = sollya_lib_cmp_in(a, b);
  sollya_lib_printf("%b in %b returns %b\n", a, b, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(res);

  sollya_lib_uninstall_msg_callback();
  sollya_lib_close();
  return 0;
}
