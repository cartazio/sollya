#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_struct_t {
  char text[32];
  int  counter;
} data_t;

int successor_impl(int a) {
  return a + 1;
}

int successor(int *res, void **args, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;

  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;

  *res = successor_impl(*((int *) (args[0])));
  return 1;
}

int successor_bis(int *res, void **args, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;

  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;

  *res = successor_impl(*((int *) (args[0])));
  return 1;
}

int strange_impl(char *str, int i, mpfr_t v) {
  sollya_lib_printf(">>>%s<<<>>>%d<<<>>>%v<<<\n", str, i, v);
  return (i == 17);
}

int strange(int *res, void **args, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;

  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;

  *res = strange_impl(((char *) (args[0])), *((int *) (args[1])), *((mpfr_t *) (args[2])));
  return 1;
}

int strange_bis(int *res, void **args, void *ptr) {
  data_t *data;

  data = (data_t *) ptr;

  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;

  *res = strange_impl(((char *) (args[0])), *((int *) (args[1])), *((mpfr_t *) (args[2])));
  return 1;
}

int empty_func(void *ptr) {
  data_t *data;

  data = (data_t *) ptr;

  sollya_lib_printf(">>>>%s<<<<>>>>%d<<<<\n", data->text, data->counter);
  (data->counter)++;

  sollya_lib_printf("External procedure only doing a side effect\n");
  return 1;
}

void dealloc(void *ptr) {
  data_t *data;
  data = (data_t *) ptr;
  sollya_lib_printf("Deallocation function called for the data pointer (%s<<<>>>%d)\n", data->text, data->counter);
  return;
}

int main(void) {
  sollya_obj_t f[19];
  sollya_externalprocedure_type_t argTypes[3];
  int i;
  char str1[1024];
  char str2[1024];
  char str3[1024];
  data_t data = { "Hello world", 0 };
  
  sollya_lib_init();

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[0] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_INTEGER, argTypes, 1, "succ", successor, &data, NULL);
  sollya_lib_printf("%b\n", f[0]);
  sollya_lib_autoprint(f[0], NULL);
  f[1] = SOLLYA_CONST_SI64(16);
  f[2] = sollya_lib_apply(f[0], f[1], NULL);
  sollya_lib_printf("%b\n", f[2]);

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[3] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_INTEGER, argTypes, 1, NULL, successor_bis, &data, NULL);
  sollya_lib_sprintf(str1, "%b", f[3]);
  sollya_lib_sprintf(str2, "proc_%p_%p", successor_bis, &data);
  sollya_lib_sprintf(str3, "%s_%p", "successor_bis", &data);
  if ((strcmp(str1, str2) == 0) || (strcmp(str1, str3) == 0)) {
    sollya_lib_printf("The behavior when no name is suggested is conform to the semantic\n");
  } else {
    sollya_lib_printf("FAILURE: the external procedure prints as \"%s\"\n",str1);
  }
  f[4] = SOLLYA_CONST_SI64(16);
  f[5] = sollya_lib_apply(f[3], f[4], NULL);
  sollya_lib_printf("%b\n", f[5]);

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_STRING;
  argTypes[1] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  argTypes[2] = SOLLYA_EXTERNALPROC_TYPE_CONSTANT;
  f[6] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_BOOLEAN, argTypes, 3, "strange_proc", strange, &data, NULL);
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
  f[11] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_BOOLEAN, argTypes, 3, NULL, strange_bis, &data, NULL);
  sollya_lib_sprintf(str1, "%b", f[11]);
  sollya_lib_sprintf(str2, "proc_%p_%p", strange_bis, &data);
  sollya_lib_sprintf(str3, "%s_%p", "strange_bis", &data);
  if ((strcmp(str1, str2) == 0) || (strcmp(str1, str3) == 0)) {
    sollya_lib_printf("The behavior when no name is suggested is conform to the semantic\n");
  } else {
    sollya_lib_printf("FAILURE: the external procedure prints as \"%s\"\n",str1);
  }
  f[12] = sollya_lib_string("Coucou");
  f[13] = SOLLYA_CONST_SI64(42);
  f[14] = SOLLYA_EXP(SOLLYA_CONST_SI64(1));
  f[15] = sollya_lib_apply(f[11], f[12], f[13], f[14], NULL);
  sollya_lib_printf("%b\n", f[15]);

  f[16] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_VOID, NULL, 0, NULL, empty_func, &data, dealloc);
  sollya_lib_sprintf(str1, "%b", f[16]);
  sollya_lib_sprintf(str2, "proc_%p_%p", empty_func, &data);
  sollya_lib_sprintf(str3, "%s_%p", "empty_func", &data);
  if ((strcmp(str1, str2) == 0) || (strcmp(str1, str3) == 0)) {
    sollya_lib_printf("The behavior when no name is suggested is conform to the semantic\n");
  } else {
    sollya_lib_printf("FAILURE: the external procedure prints as \"%s\"\n",str1);
  }
  f[17] = sollya_lib_apply(f[16], NULL);
  f[18] = sollya_lib_copy_obj(f[16]);
  sollya_lib_clear_obj(f[16]);
  sollya_lib_printf("Deallocation function must not yet have been called\n");
  f[16] = SOLLYA_X_;

  for(i=0;i<=18;i++) sollya_lib_clear_obj(f[i]);

  sollya_lib_close();

  return 0;
}
