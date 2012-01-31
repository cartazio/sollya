#include <sollya.h>

int main(void) {
  sollya_obj_t a[6], b[6];
  int i;

  sollya_lib_init();

  a[0] = SOLLYA_POW(SOLLYA_ADD(SOLLYA_POW(SOLLYA_ADD(SOLLYA_MUL(SOLLYA_CONST(2.0),SOLLYA_X_),SOLLYA_CONST(3.0)),SOLLYA_CONST(2.0)),SOLLYA_POW(SOLLYA_X_,SOLLYA_CONST(4.0))),SOLLYA_CONST(2.0));
  a[1] = SOLLYA_EXP(SOLLYA_X_);
  a[2] = SOLLYA_ADD(sollya_lib_copy_obj(a[0]),SOLLYA_SIN(sollya_lib_copy_obj(a[0])));
  a[3] = sollya_lib_parse_string("horner((x^4 + x^2 + ((x + 3)^2)/(x + 4))^7);");
  a[4] = sollya_lib_parse_string("(((x + 1)/(x + 2)) * ((x + 3)/(x + 4)))/(((x + 5)/(x + 6)) * ((x + 7)/(x + 8)));");
  a[5] = sollya_lib_parse_string("(((x + 1)/(x + 2)) * ((x + 3)/(x + 4)))/x;");
    
  for (i=0;i<6;i++) {
    b[i] = sollya_lib_expand(a[i]);
  }

  for (i=0;i<6;i++) {
    sollya_lib_printf("%b with all polynomial subtrees maximally expanded is %b.\n",a[i],b[i]);
  }

  for (i=0;i<6;i++) {
    sollya_lib_clear_obj(a[i]);
    sollya_lib_clear_obj(b[i]);
  }
  
  sollya_lib_close();
  return 0;
}
