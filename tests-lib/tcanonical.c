#include <sollya.h>

int main(void) {
  sollya_obj_t a[4], b[4];
  int i;

  sollya_lib_init();

  a[0] = SOLLYA_POW(SOLLYA_ADD(SOLLYA_POW(SOLLYA_ADD(SOLLYA_MUL(SOLLYA_CONST(2.0),SOLLYA_X_),SOLLYA_CONST(3.0)),SOLLYA_CONST(2.0)),SOLLYA_POW(SOLLYA_X_,SOLLYA_CONST(4.0))),SOLLYA_CONST(2.0));
  a[1] = SOLLYA_EXP(SOLLYA_X_);
  a[2] = SOLLYA_ADD(sollya_lib_copy_obj(a[0]),SOLLYA_SIN(sollya_lib_copy_obj(a[0])));
  a[3] = sollya_lib_parse_string("horner((x^4 + x^2 + 1)^23);");
    
  for (i=0;i<4;i++) {
    b[i] = sollya_lib_canonical(a[i]);
  }

  for (i=0;i<4;i++) {
    sollya_lib_printf("%b with all polynomial subtrees written in the canonical basis is %b.\n",a[i],b[i]);
  }

  for (i=0;i<4;i++) {
    sollya_lib_clear_obj(a[i]);
    sollya_lib_clear_obj(b[i]);
  }
  
  sollya_lib_close();
  return 0;
}
