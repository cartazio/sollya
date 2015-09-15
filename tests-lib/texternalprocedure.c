#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <string.h>

int successor_impl(int a) {
  return a + 1;
}

int successor(int *res, void **args) {
  *res = successor_impl(*((int *) (args[0])));
  return 1;
}

int successor_bis(int *res, void **args) {
  *res = successor_impl(*((int *) (args[0])));
  return 1;
}

int strange_impl(char *str, int i, mpfr_t v) {
  sollya_lib_printf(">>>%s<<<>>>%d<<<>>>%v<<<\n", str, i, v);
  return (i == 17);
}

int strange(int *res, void **args) {
  *res = strange_impl(((char *) (args[0])), *((int *) (args[1])), *((mpfr_t *) (args[2])));
  return 1;
}

int strange_bis(int *res, void **args) {
  *res = strange_impl(((char *) (args[0])), *((int *) (args[1])), *((mpfr_t *) (args[2])));
  return 1;
}

int main(void) {
  sollya_obj_t f[11];
  sollya_externalprocedure_type_t argTypes[3];
  int i;

  sollya_lib_init();

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[0] = sollya_lib_externalprocedure(SOLLYA_EXTERNALPROC_TYPE_INTEGER, argTypes, 1, "succ", successor);
  sollya_lib_printf("%b\n", f[0]);
  sollya_lib_autoprint(f[0], NULL);
  f[1] = SOLLYA_CONST_SI64(16);
  f[2] = sollya_lib_apply(f[0], f[1], NULL);
  sollya_lib_printf("%b\n", f[2]);

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[3] = sollya_lib_externalprocedure(SOLLYA_EXTERNALPROC_TYPE_INTEGER, argTypes, 1, NULL, successor_bis);
  sollya_lib_printf("%b\n", f[3]);
  sollya_lib_autoprint(f[3], NULL);
  f[4] = SOLLYA_CONST_SI64(16);
  f[5] = sollya_lib_apply(f[3], f[4], NULL);
  sollya_lib_printf("%b\n", f[5]);

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_STRING;
  argTypes[1] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  argTypes[2] = SOLLYA_EXTERNALPROC_TYPE_CONSTANT;
  f[6] = sollya_lib_externalprocedure(SOLLYA_EXTERNALPROC_TYPE_BOOLEAN, argTypes, 3, "strange_proc", strange);
  sollya_lib_printf("%b\n", f[6]);
  sollya_lib_autoprint(f[6], NULL);
  f[7] = sollya_lib_string("Coucou");
  f[8] = SOLLYA_CONST_SI64(17);
  f[9] = SOLLYA_EXP(SOLLYA_CONST_SI64(1));
  f[10] = sollya_lib_apply(f[6], f[7], f[8], f[9], NULL);
  sollya_lib_printf("%b\n", f[10]);

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_STRING;
  argTypes[1] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  argTypes[2] = SOLLYA_EXTERNALPROC_TYPE_CONSTANT;
  f[6] = sollya_lib_externalprocedure(SOLLYA_EXTERNALPROC_TYPE_BOOLEAN, argTypes, 3, NULL, strange_bis);
  sollya_lib_printf("%b\n", f[6]);
  sollya_lib_autoprint(f[6], NULL);
  f[7] = sollya_lib_string("Coucou");
  f[8] = SOLLYA_CONST_SI64(42);
  f[9] = SOLLYA_EXP(SOLLYA_CONST_SI64(1));
  f[10] = sollya_lib_apply(f[6], f[7], f[8], f[9], NULL);
  sollya_lib_printf("%b\n", f[10]);
  
  for(i=0;i<=10;i++) sollya_lib_clear_obj(f[i]);
  
  sollya_lib_close();
  
  return 0;
}
