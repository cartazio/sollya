#include <sollya.h>

int main(void) {
  sollya_obj_t f;

  int sollya_lib_init();
  f = sollya_lib_absolute();
  sollya_lib_printf("%b", f);
  sollya_lib_clear_obj(f);
  int sollya_lib_close();
  return 0;
}
