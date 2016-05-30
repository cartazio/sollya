#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_struct_t {
  char text[32];
  int  counter;
} data_t;

char *externalprocTypeToString(sollya_externalprocedure_type_t t) {
  switch (t) {
  case SOLLYA_EXTERNALPROC_TYPE_VOID:
    return "void";
  case SOLLYA_EXTERNALPROC_TYPE_CONSTANT:
    return "constant";
  case SOLLYA_EXTERNALPROC_TYPE_FUNCTION:
    return "function";
  case SOLLYA_EXTERNALPROC_TYPE_RANGE:
    return "range";
  case SOLLYA_EXTERNALPROC_TYPE_INTEGER:
    return "integer";
  case SOLLYA_EXTERNALPROC_TYPE_STRING:
    return "string";
  case SOLLYA_EXTERNALPROC_TYPE_BOOLEAN:
    return "boolean";
  case SOLLYA_EXTERNALPROC_TYPE_OBJECT:
    return "object";
  case SOLLYA_EXTERNALPROC_TYPE_CONSTANT_LIST:
    return "list of constant";
  case SOLLYA_EXTERNALPROC_TYPE_FUNCTION_LIST:
    return "list of function";
  case SOLLYA_EXTERNALPROC_TYPE_RANGE_LIST:
    return "list of range";
  case SOLLYA_EXTERNALPROC_TYPE_INTEGER_LIST:
    return "list of integer";
  case SOLLYA_EXTERNALPROC_TYPE_STRING_LIST:
    return "list of string";
  case SOLLYA_EXTERNALPROC_TYPE_BOOLEAN_LIST:
    return "list of boolean";
  case SOLLYA_EXTERNALPROC_TYPE_OBJECT_LIST:
    return "list of object";
  default:
    break;
  }
  return "unknown";
}

void dealloc(void *data) {
  (void) data;
  return;
}

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

int empty_func2() {
  sollya_lib_printf("External procedure only doing a side effect number 2\n");
  return 1;
}

void dealloc2(void *ptr) {
  data_t *data;
  data = (data_t *) ptr;
  sollya_lib_printf("Deallocation function called for the data pointer (%s<<<>>>%d)\n", data->text, data->counter);
  return;
}

int main(void) {
  sollya_obj_t f[20];
  sollya_externalprocedure_type_t argTypes[3];
  sollya_externalprocedure_type_t result_type;
  sollya_externalprocedure_type_t *argument_types;
  int res;
  int arity;
  int i;
  char str1[1024];
  char str2[1024];
  char str3[1024];
  data_t data = { "Hello world", 0 };
  data_t tryData = { "Grias Eahna Good.", 17 };
  void *func;
  void *resData;
  int myInt;
  int result;
  void *args[3];
  mpfr_t myMpfr;
  void (*resDealloc)(void *);
  
  sollya_lib_init();

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[0] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_INTEGER, argTypes, 1, "succ", successor, &data, NULL);
  sollya_lib_printf("%b\n", f[0]);
  sollya_lib_autoprint(f[0], NULL);
  f[1] = SOLLYA_CONST_SI64(16);
  f[2] = sollya_lib_apply(f[0], f[1], NULL);
  sollya_lib_printf("%b\n", f[2]);

  result_type = (sollya_externalprocedure_type_t) -1; argument_types = NULL; arity = -1; func = NULL; resData = NULL; resDealloc = NULL;
  res = sollya_lib_decompose_externalprocedure_with_data(&result_type, &argument_types, &arity, &func, &resData, &resDealloc, f[0]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data has worked on \"%b\": arity = %d, result type = %s, argument types = ", f[0], arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    sollya_lib_printf("The function pointer is %s\n", ((((void *) func) == ((void *) successor)) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((resData == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == NULL) ? "okay" : "wrong"));
    if ((arity == 1) &&
	(result_type == SOLLYA_EXTERNALPROC_TYPE_INTEGER) &&
	(argument_types[0] == SOLLYA_EXTERNALPROC_TYPE_INTEGER)) {
      sollya_lib_printf("Trying out the function\n");
      myInt = 42;
      args[0] = &myInt;
      res = ((int (*)(int *, void **, void *)) func)(&result, args, &tryData);
      sollya_lib_printf("The function has signaled %s, the result is %d, the data counter is %d\n", (res ? "success" : "failure"), result, tryData.counter);
    }    
    sollya_lib_free(argument_types);
  }

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[3] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_INTEGER, argTypes, 1, NULL, successor_bis, &data, dealloc);
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

  result_type = (sollya_externalprocedure_type_t) -1; argument_types = NULL; arity = -1; func = NULL; resData = NULL; resDealloc = NULL;
  res = sollya_lib_decompose_externalprocedure_with_data(&result_type, &argument_types, &arity, &func, &resData, &resDealloc, f[3]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data has worked: arity = %d, result type = %s, argument types = ", arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    sollya_lib_printf("The function pointer is %s\n", ((((void *) func) == ((void *) successor_bis)) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((resData == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
    if ((arity == 1) &&
	(result_type == SOLLYA_EXTERNALPROC_TYPE_INTEGER) &&
	(argument_types[0] == SOLLYA_EXTERNALPROC_TYPE_INTEGER)) {
      sollya_lib_printf("Trying out the function\n");
      myInt = 42;
      args[0] = &myInt;
      res = ((int (*)(int *, void **, void *)) func)(&result, args, &tryData);
      sollya_lib_printf("The function has signaled %s, the result is %d, the data counter is %d\n", (res ? "success" : "failure"), result, tryData.counter);
    }    
    sollya_lib_free(argument_types);
  }

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_STRING;
  argTypes[1] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  argTypes[2] = SOLLYA_EXTERNALPROC_TYPE_CONSTANT;
  f[6] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_BOOLEAN, argTypes, 3, "strange_proc", strange, &data, dealloc);
  sollya_lib_printf("%b\n", f[6]);
  sollya_lib_autoprint(f[6], NULL);
  f[7] = sollya_lib_string("Coucou");
  f[8] = SOLLYA_CONST_SI64(17);
  f[9] = SOLLYA_EXP(SOLLYA_CONST_SI64(1));
  f[10] = sollya_lib_apply(f[6], f[7], f[8], f[9], NULL);
  sollya_lib_printf("%b\n", f[10]);

  result_type = (sollya_externalprocedure_type_t) -1; argument_types = NULL; arity = -1; func = NULL; resData = NULL; resDealloc = NULL;
  res = sollya_lib_decompose_externalprocedure_with_data(&result_type, &argument_types, &arity, &func, &resData, &resDealloc, f[6]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data has worked on \"%b\": arity = %d, result type = %s, argument types = ", f[6], arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    sollya_lib_printf("The function pointer is %s\n", ((((void *) func) == ((void *) strange)) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((resData == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
    if ((arity == 3) &&
	(result_type == SOLLYA_EXTERNALPROC_TYPE_BOOLEAN) &&
	(argument_types[0] == SOLLYA_EXTERNALPROC_TYPE_STRING) &&
	(argument_types[1] == SOLLYA_EXTERNALPROC_TYPE_INTEGER) &&
	(argument_types[2] == SOLLYA_EXTERNALPROC_TYPE_CONSTANT)) {
      sollya_lib_printf("Trying out the function\n");
      args[0] = (void *) "Hello world";
      myInt = 42;
      args[1] = &myInt;
      mpfr_init2(myMpfr, 53);
      mpfr_set_si(myMpfr, 1664, GMP_RNDN);
      args[2] = &myMpfr;
      res = ((int (*)(int *, void **, void *)) func)(&result, args, &tryData);
      sollya_lib_printf("The function has signaled %s, the result is %s, the data counter is %d\n", (res ? "success" : "failure"), (result ? "true" : "false"), tryData.counter);
      mpfr_clear(myMpfr);
    }    
    sollya_lib_free(argument_types);
  }

  argTypes[0] = SOLLYA_EXTERNALPROC_TYPE_STRING;
  argTypes[1] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  argTypes[2] = SOLLYA_EXTERNALPROC_TYPE_CONSTANT;
  f[11] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_BOOLEAN, argTypes, 3, NULL, strange_bis, &data, dealloc);
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

  result_type = (sollya_externalprocedure_type_t) -1; argument_types = NULL; arity = -1; func = NULL; resData = NULL; resDealloc = NULL;
  res = sollya_lib_decompose_externalprocedure_with_data(&result_type, &argument_types, &arity, &func, &resData, &resDealloc, f[11]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data has worked: arity = %d, result type = %s, argument types = ", arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    sollya_lib_printf("The function pointer is %s\n", ((((void *) func) == ((void *) strange_bis)) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((resData == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc)) ? "okay" : "wrong"));
    if ((arity == 3) &&
	(result_type == SOLLYA_EXTERNALPROC_TYPE_BOOLEAN) &&
	(argument_types[0] == SOLLYA_EXTERNALPROC_TYPE_STRING) &&
	(argument_types[1] == SOLLYA_EXTERNALPROC_TYPE_INTEGER) &&
	(argument_types[2] == SOLLYA_EXTERNALPROC_TYPE_CONSTANT)) {
      sollya_lib_printf("Trying out the function\n");
      args[0] = (void *) "Hello world";
      myInt = 42;
      args[1] = &myInt;
      mpfr_init2(myMpfr, 53);
      mpfr_set_si(myMpfr, 1664, GMP_RNDN);
      args[2] = &myMpfr;
      res = ((int (*)(int *, void **, void *)) func)(&result, args, &tryData);
      sollya_lib_printf("The function has signaled %s, the result is %s, the data counter is %d\n", (res ? "success" : "failure"), (result ? "true" : "false"), tryData.counter);
      mpfr_clear(myMpfr);
    }    
    sollya_lib_free(argument_types);
  }

  f[16] = sollya_lib_externalprocedure_with_data(SOLLYA_EXTERNALPROC_TYPE_VOID, NULL, 0, NULL, empty_func, &data, dealloc2);
  f[17] = sollya_lib_apply(f[16], NULL);
  sollya_lib_printf("%b\n", f[17]);

  result_type = (sollya_externalprocedure_type_t) -1; argument_types = NULL; arity = -1; func = NULL; resData = NULL; resDealloc = NULL;
  res = sollya_lib_decompose_externalprocedure_with_data(&result_type, &argument_types, &arity, &func, &resData, &resDealloc, f[16]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data has worked: arity = %d, result type = %s, argument types = ", arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    sollya_lib_printf("The function pointer is %s\n", ((((void *) func) == ((void *) empty_func)) ? "okay" : "wrong"));
    sollya_lib_printf("The data pointer is %s\n", ((resData == ((void *) (&data))) ? "okay" : "wrong"));
    sollya_lib_printf("The deallocation function pointer is %s\n", ((((void *) resDealloc) == ((void *) dealloc2)) ? "okay" : "wrong"));
  }

  f[18] = sollya_lib_externalprocedure(SOLLYA_EXTERNALPROC_TYPE_VOID, NULL, 0, NULL, empty_func2);
  f[19] = sollya_lib_apply(f[18], NULL);
  sollya_lib_printf("%b\n", f[19]);

  result_type = (sollya_externalprocedure_type_t) -1; argument_types = NULL; arity = -1; func = NULL; resData = NULL; resDealloc = NULL;
  res = sollya_lib_decompose_externalprocedure_with_data(&result_type, &argument_types, &arity, &func, &resData, &resDealloc, f[18]);
  if (!res) sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data did not work (as expected) on something not constructed with sollya_lib_externalprocedure_with_data\n");
  else sollya_lib_printf("sollya_lib_decompose_externalprocedure_with_data worked in an unexpected case\n"); 

  for(i=0;i<=19;i++) sollya_lib_clear_obj(f[i]);
  
  sollya_lib_close();
  
  return 0;
}
