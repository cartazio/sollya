#include <sollya.h>

int callback(int message) {
  if (message == SOLLYA_MSG_EXPR_NOT_CORRECTLY_TYPED)
    sollya_lib_printf("The following test produces a typing error\n");
  else
    sollya_lib_printf("Unexpected message\n");
  return 0;

}
int main(void) {
  sollya_obj_t true, false, a, res;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback);

  true = sollya_lib_true();
  false = sollya_lib_false();

  a = sollya_lib_copy_obj(true);
  res = sollya_lib_negate(a);
  sollya_lib_printf("the negation of %b gives %b\n", a, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(res);

  a = sollya_lib_copy_obj(false);
  res = sollya_lib_negate(a);
  sollya_lib_printf("the negation of %b gives %b\n", a, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(res);

  a = SOLLYA_CONST(17);
  res = sollya_lib_negate(a);
  sollya_lib_printf("the negation of %b gives %b\n", a, res);
  sollya_lib_clear_obj(a);
  sollya_lib_clear_obj(res);


  sollya_lib_clear_obj(true);
  sollya_lib_clear_obj(false);
  sollya_lib_uninstall_msg_callback();
  sollya_lib_close();
  return 0;
}