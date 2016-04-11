#include <sollya.h>


sollya_obj_t stupid_wrapper(sollya_obj_t arg1, ...) {
  va_list va;
  sollya_obj_t a;
  va_start(va, arg1);
  a = sollya_lib_v_execute_procedure(arg1, va);
  va_end(va);
  return a;
}

int callback(sollya_msg_t msg, void *data) {
  (void)msg; /* Avoiding "unused parameter" warning */
  (void)data; /* Avoiding "unused parameter" warning */
  return 0;
}

int main(void) {
  sollya_obj_t proc1, proc2, proc3, one, two, three, proc1Res, proc1Resb, proc1Resc, proc2Res, proc3Res, proc3Resb;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback, NULL);

  proc1 = sollya_lib_parse_string("proc (a,b) { \"a = \", a, \", b = \", b; return a + b; }");
  proc2 = sollya_lib_parse_string("proc () { write(\"Coucou\\n\"); }");
  proc3 = sollya_lib_parse_string("proc (l = ...) { \"l = \", l; return exp(_x_); }");

  one = SOLLYA_CONST_SI64(1);
  two = SOLLYA_CONST_SI64(2);
  three = SOLLYA_CONST_SI64(3);

  proc1Res = stupid_wrapper(proc1, one, two, NULL);
  proc1Resb = stupid_wrapper(proc1, one, NULL);
  proc1Resc = stupid_wrapper(proc1, one, two, three, NULL);
  proc2Res = stupid_wrapper(proc2, NULL);
  proc3Res = stupid_wrapper(proc3, one, two, three, NULL);
  proc3Resb = stupid_wrapper(proc3, NULL);

  sollya_lib_printf("proc1 = %b, args = %b, %b: result = %b\n", proc1, one, two, proc1Res);
  sollya_lib_printf("proc1 = %b, args = %b: result = %b\n", proc1, one, proc1Resb);
  sollya_lib_printf("proc1 = %b, args = %b, %b, %b: result = %b\n", proc1, one, two, three, proc1Resc);
  sollya_lib_printf("proc2 = %b: result = %b\n", proc2, proc2Res);
  sollya_lib_printf("proc3 = %b, args = %b, %b, %b: result = %b\n", proc3, one, two, three, proc3Res);
  sollya_lib_printf("proc3 = %b: result = %b\n", proc3, one, two, three, proc3Resb);

  sollya_lib_clear_obj(proc3Res);
  sollya_lib_clear_obj(proc3Resb);
  sollya_lib_clear_obj(proc2Res);
  sollya_lib_clear_obj(proc1Res);
  sollya_lib_clear_obj(proc1Resb);
  sollya_lib_clear_obj(proc1Resc);
  sollya_lib_clear_obj(three);
  sollya_lib_clear_obj(two);
  sollya_lib_clear_obj(one);
  sollya_lib_clear_obj(proc3);
  sollya_lib_clear_obj(proc2);
  sollya_lib_clear_obj(proc1);

  sollya_lib_close();
  return 0;
}
