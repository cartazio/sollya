#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>

void euler_gamma(mpfr_t res, mp_prec_t prec) {
  mpfr_set_prec(res, prec);
  mpfr_const_euler(res, GMP_RNDN);
  return;
}

int stupid1(mpfi_t result, mpfi_t x, int n) {
  (void)x; /* Avoiding "unused parameter" warning */
  (void)n; /* Avoiding "unused parameter" warning */
  mpfi_set_ui(result, 0);
  return 0;
}
char *type_to_string(sollya_base_function_t type) {
  switch (type) {
  case SOLLYA_BASE_FUNC_ASINH: return "SOLLYA_BASE_FUNC_ASINH";
  case SOLLYA_BASE_FUNC_LIBRARYCONSTANT: return "SOLLYA_BASE_FUNC_LIBRARYCONSTANT";
  case SOLLYA_BASE_FUNC_CONSTANT: return "SOLLYA_BASE_FUNC_CONSTANT";
  case SOLLYA_BASE_FUNC_PI: return "SOLLYA_BASE_FUNC_PI";
  case SOLLYA_BASE_FUNC_LIBRARYFUNCTION: return "SOLLYA_BASE_FUNC_LIBRARYFUNCTION";
  case SOLLYA_BASE_FUNC_PROCEDUREFUNCTION: return "SOLLYA_BASE_FUNC_PROCEDUREFUNCTION";
  case SOLLYA_BASE_FUNC_POW: return "SOLLYA_BASE_FUNC_POW";
  default: return "UNKNOWN";
  }
}

int main(void) {
  sollya_base_function_t type;
  sollya_obj_t f, obj1, obj2, obj3;
  int res;

  sollya_lib_init();

  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_ASINH;
  obj1 = SOLLYA_X_;
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_ASINH;
  obj1 = SOLLYA_X_; obj2 = sollya_lib_parse_string("[|1,2|]");
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_ASINH;
  obj1 = sollya_lib_parse_string("[1,2]"); obj2 = SOLLYA_X_;
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_ASINH;
  obj1 = sollya_lib_parse_string("[|1,2|]"); obj2 = SOLLYA_X_;
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_LIBRARYCONSTANT;
  obj1 = SOLLYA_X_;
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_LIBRARYCONSTANT;
  obj1 = sollya_lib_parse_string("[|1,2|]");
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_LIBRARYCONSTANT;
  obj1 = sollya_lib_libraryconstant("superconst", euler_gamma);
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);

  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_CONSTANT;
  obj1 = sollya_lib_libraryconstant("superconst", euler_gamma);
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_CONSTANT;
  obj1 = SOLLYA_CONST(3.0);
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_PI;
  res = sollya_lib_construct_function(&f, type, NULL);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_PI;
  obj1 = SOLLYA_CONST(3.0);
  res = sollya_lib_construct_function(&f, type, obj1);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_LIBRARYFUNCTION;
  obj1 = SOLLYA_CONST(3.0); obj2 = sollya_lib_parse_string("[|1,2|]");
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_LIBRARYFUNCTION;
  obj1 = SOLLYA_CONST(3.0); obj2 = sollya_lib_build_function_libraryfunction(SOLLYA_X_, "superfunc", stupid1);
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_LIBRARYFUNCTION;
  obj1 = SOLLYA_CONST(3.0); obj2 = sollya_lib_build_function_libraryfunction(SOLLYA_COS(SOLLYA_X_), "superfunc", stupid1);
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_PROCEDUREFUNCTION;
  obj1 = sollya_lib_parse_string("[1,2]"); obj2 = sollya_lib_build_function_libraryfunction(SOLLYA_COS(SOLLYA_X_), "superfunc", stupid1);
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_PROCEDUREFUNCTION;
  obj1 = sollya_lib_parse_string("[1,2]");
  obj2 = sollya_lib_parse_string("proc(X,n,p) {var res, oldPrec; oldPrec = prec; prec = p!; res = exp(X); prec = oldPrec!; return res; };");
  obj2 = sollya_lib_build_function_procedurefunction(SOLLYA_X_, obj2);
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_PROCEDUREFUNCTION;
  obj1 = SOLLYA_CONST(3.0);
  obj2 = sollya_lib_parse_string("proc(X,n,p) {var res, oldPrec; oldPrec = prec; prec = p!; res = exp(X); prec = oldPrec!; return res; };");
  obj2 = sollya_lib_build_function_procedurefunction(SOLLYA_COS(SOLLYA_X_), obj2);
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);


  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_POW;
  obj1 = SOLLYA_X_; obj2 = SOLLYA_CONST(3.0);
  res = sollya_lib_construct_function(&f, type, obj1, obj2);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);

  f = NULL; obj1 = NULL; obj2 = NULL; obj3 = NULL; res = -1;
  type = SOLLYA_BASE_FUNC_POW;
  obj1 = SOLLYA_X_; obj2 = SOLLYA_X_; obj3 = sollya_lib_parse_string("[|1,2|]");
  res = sollya_lib_construct_function(&f, type, obj1, obj2, obj3);
  if (res)
    sollya_lib_printf("Constructed %b from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  else
    sollya_lib_printf("Impossible to construct (f=%b) from type=%s, obj1=%b, obj2=%b and obj3=%b\n", f, type_to_string(type), obj1, obj2, obj3);
  if (f != NULL) sollya_lib_clear_obj(f);
  if (obj1 != NULL) sollya_lib_clear_obj(obj1);
  if (obj2 != NULL) sollya_lib_clear_obj(obj2);
  if (obj3 != NULL) sollya_lib_clear_obj(obj3);



  sollya_lib_close();
  return 0;
}
