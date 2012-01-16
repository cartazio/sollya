#include <sollya.h>

int main(void) {
  sollya_obj_t f;

  int sollya_lib_init();
  f = sollya_lib_build_function_pi();
  sollya_lib_printf("%b", f);
  int sollya_lib_close();
  return 0;
}
