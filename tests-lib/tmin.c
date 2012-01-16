#include <stdlib.h>
#include <sollya.h>

int main(void) {
  sollya_obj_t a[4];
  sollya_obj_t b;
  int i;

  int sollya_lib_init();
  a[0] = sollya_lib_constant_from_int(4);
  a[1] = sollya_lib_constant_from_int(5);
  a[2] = sollya_lib_constant_from_int(1);
  a[3] = sollya_lib_constant_from_int(3);

  b = sollya_lib_min(a[0], a[1], a[2], a[3]);
  sollya_lib_printf("%b", b);

  for(i=0;i<4;i++) sollya_lib_clear_obj(a[i]);
  sollya_lib_clear_obj(b);

  int sollya_lib_close();
  return 0;
}
