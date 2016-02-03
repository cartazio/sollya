#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>

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

int succ_impl(int a) {
  return a + 1;
}

int succ(int *res, void **args) {
  *res = succ_impl(*((int *) (args[0])));
  return 1;
}

int main(void) {
  int i;
  sollya_obj_t f[4];
  int res;
  sollya_externalprocedure_type_t result_type;
  sollya_externalprocedure_type_t *argument_types;
  sollya_externalprocedure_type_t myArgTypes[1];
  int arity;
  void *func;
  void *args[3];
  int result;
  int myInt;
  mpfr_t myMpfr;

  sollya_lib_init();

  f[0] = SOLLYA_COS(SOLLYA_X_);
  f[1] = sollya_lib_parse_string("(proc() { externalproc(funny, \"./libraryexample.a\", (string, integer, constant) -> boolean); return funny;})()");
  f[2] = sollya_lib_parse_string("(proc() { externalproc(zeit, \"./libraryexample.a\", void -> integer); return zeit;})()");
  myArgTypes[0] = SOLLYA_EXTERNALPROC_TYPE_INTEGER;
  f[3] = sollya_lib_externalprocedure(SOLLYA_EXTERNALPROC_TYPE_INTEGER, myArgTypes, 1, "succ", succ);

  argument_types = NULL;
  res = sollya_lib_decompose_externalprocedure(&result_type, &argument_types, &arity, &func, f[0]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure has signaled success on \"%b\" while it should not have worked.\n", f[0]);
    sollya_lib_free(argument_types);
  } else {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure correctly signals that \"%b\" cannot be decomposed as an external function.\n", f[0]);
  }

  argument_types = NULL;
  res = sollya_lib_decompose_externalprocedure(&result_type, &argument_types, &arity, &func, f[1]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure has worked on \"%b\": arity = %d, result type = %s, argument types = ", f[1], arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
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
      res = ((int (*)(int *, void **)) func)(&result, args);
      sollya_lib_printf("The function has signaled %s, the result is %s\n", (res ? "success" : "failure"), (result ? "true" : "false"));
      mpfr_clear(myMpfr);
    }    
    sollya_lib_free(argument_types);
  }

  argument_types = NULL;
  res = sollya_lib_decompose_externalprocedure(&result_type, &argument_types, &arity, &func, f[2]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure has worked on \"%b\": arity = %d, result type = %s, argument types = ", f[2], arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    if ((arity == 0) &&
	(result_type == SOLLYA_EXTERNALPROC_TYPE_INTEGER)) {
      sollya_lib_printf("Trying out the function\n");
      result = 17;
      res = ((int (*)(int *)) func)(&result);
      sollya_lib_printf("The function has signaled %s, the result is %s\n", (res ? "success" : "failure"), ((result != 17) ? "okay" : "not okay"));
    }    
    sollya_lib_free(argument_types);
  }

  argument_types = NULL;
  res = sollya_lib_decompose_externalprocedure(&result_type, &argument_types, &arity, &func, f[3]);
  if (res) {
    sollya_lib_printf("sollya_lib_decompose_externalprocedure has worked on \"%b\": arity = %d, result type = %s, argument types = ", f[3], arity, externalprocTypeToString(result_type));
    for (i=0;i<arity;i++) {
      sollya_lib_printf("%s ", externalprocTypeToString(argument_types[i]));
    }
    sollya_lib_printf("\n");
    if ((arity == 1) &&
	(result_type == SOLLYA_EXTERNALPROC_TYPE_INTEGER) &&
	(argument_types[0] == SOLLYA_EXTERNALPROC_TYPE_INTEGER)) {
      sollya_lib_printf("Trying out the function\n");
      myInt = 41;
      args[0] = &myInt;
      res = ((int (*)(int *, void **)) func)(&result, args);
      sollya_lib_printf("The function has signaled %s, the result is %d\n", (res ? "success" : "failure"), result);
    }    
    sollya_lib_free(argument_types);
  }
  
  for(i=0;i<=3;i++)  sollya_lib_clear_obj(f[i]);

  sollya_lib_close();
  return 0;
}
