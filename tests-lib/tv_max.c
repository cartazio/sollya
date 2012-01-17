#include <sollya.h>
#include <sollya.h>
#include <stdarg.h>

sollya_obj_t stupid_wrapper(sollya_obj_t arg1, ...) {
  va_list va;
  sollya_obj_t a;
  va_start(va, arg1);
  a = sollya_lib_v_max(arg1, va);
  va_end(va);
  return a;
}

int main(void) {
  sollya_obj_t a[4];
  sollya_obj_t b,c;
  int i;

  sollya_lib_init();

  /* Tests a simple maximum */
  a[0] = sollya_lib_constant_from_int(4);
  a[1] = sollya_lib_constant_from_int(5);
  a[2] = sollya_lib_constant_from_int(1);
  a[3] = sollya_lib_constant_from_int(3);

  b = sollya_lib_max(a[0], a[1], a[2], a[3], NULL);
  sollya_lib_printf("max(4,5,1,3) returns %b\n", b);
  sollya_lib_clear_obj(b);

  c = sollya_lib_list(a, 4);
  b = sollya_lib_max(c, NULL);
  sollya_lib_printf("max(%b) returns %b\n", c, b);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);

  for(i=0;i<4;i++) sollya_lib_clear_obj(a[i]);


  /* Tests a tricky case where the maximum is impossible to detect */
  a[0] = sollya_lib_parse_string("17 + log2(13)/log2(9);");
  a[1] = sollya_lib_parse_string("17 + log(13)/log(9);");

  b = sollya_lib_max(a[0], a[1], NULL);
  sollya_lib_printf("max of 17 + log2(13)/log2(9) and 17 + log(13)/log(9) returns %b\n", b);
  sollya_lib_clear_obj(b);

  c = sollya_lib_list(a, 2);
  b = sollya_lib_max(c, NULL);
  sollya_lib_printf("max(%b) returns %b\n", c, b);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);

  sollya_lib_clear_obj(a[0]);
  sollya_lib_clear_obj(a[1]);


  /* Tests what happens when a NaN is in the list */
  a[0] = sollya_lib_constant_from_int(2);
  a[1] = sollya_lib_parse_string("NaN;");
  a[2] = sollya_lib_constant_from_int(1);

  b = sollya_lib_max(a[0], a[1], a[2], NULL);
  sollya_lib_printf("max(2,NaN,1) returns %b\n", b);
  sollya_lib_clear_obj(b);

  c = sollya_lib_list(a, 3);
  b = sollya_lib_max(c, NULL);
  sollya_lib_printf("max(%b) returns %b\n", c, b);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);

  sollya_lib_clear_obj(a[0]);
  sollya_lib_clear_obj(a[1]);
  sollya_lib_clear_obj(a[2]);

  a[0] = sollya_lib_parse_string("NaN;");
  a[1] = sollya_lib_constant_from_int(1);
  a[2] = sollya_lib_constant_from_int(2);

  b = sollya_lib_max(a[0], a[1], a[2], NULL);
  sollya_lib_printf("max(1,NaN,2) returns %b\n", b);
  sollya_lib_clear_obj(b);

  c = sollya_lib_list(a, 3);
  b = sollya_lib_max(c, NULL);
  sollya_lib_printf("max(%b) returns %b\n", c, b);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);

  sollya_lib_clear_obj(a[0]);
  sollya_lib_clear_obj(a[1]);
  sollya_lib_clear_obj(a[2]);

  /* Tests maximum of only one element */
  a[0] = sollya_lib_constant_from_int(17);
  b = sollya_lib_max(a[0], NULL);
  sollya_lib_printf("max of 17 returns %b\n", b);
  sollya_lib_clear_obj(b);

  c = sollya_lib_list(a, 1);
  b = sollya_lib_max(c, NULL);
  sollya_lib_printf("max(%b) returns %b\n", c, b);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);

  sollya_lib_clear_obj(a[0]);

  /* Tests maximum of an empty list */
  c = sollya_lib_list(NULL, 0);
  b = sollya_lib_max(c, NULL);
  sollya_lib_printf("max of an empty list returns %b\n", b);
  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);


  /* Tests what happens if a list is given, followed by a constant */
  a[0] = sollya_lib_constant_from_int(4);
  a[1] = sollya_lib_constant_from_int(1);
  a[2] = sollya_lib_constant_from_int(3);
  a[3] = sollya_lib_constant_from_int(5);
  c = sollya_lib_list(a, 3);
  b = sollya_lib_max(c, a[3], NULL);
  sollya_lib_printf("max(%b, 5) returns %b\n", c, b);

  sollya_lib_clear_obj(b);
  sollya_lib_clear_obj(c);
  for(i=0;i<4;i++) sollya_lib_clear_obj(a[i]);


  sollya_lib_close();
  return 0;
}