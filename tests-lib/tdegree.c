#include <sollya.h>

int main(void) {
  sollya_obj_t a[7], b[7];
  int i;

  sollya_lib_init();

  a[0] = SOLLYA_CONST(1.0);
  a[1] = SOLLYA_X_;
  a[2] = SOLLYA_POW(SOLLYA_X_,SOLLYA_CONST(2.0));
  a[3] = SOLLYA_MUL(SOLLYA_X_,SOLLYA_X_);
  a[4] = SOLLYA_POW(SOLLYA_ADD(SOLLYA_POW(SOLLYA_X_,SOLLYA_CONST(2.0)),SOLLYA_CONST(1.0)),SOLLYA_CONST(5.0));
  a[5] = SOLLYA_CONST(0.0);
  a[6] = SOLLYA_EXP(SOLLYA_X_);

  for (i=0;i<7;i++) {
    b[i] = sollya_lib_degree(a[i]);
  }

  for (i=0;i<7;i++) {
    sollya_lib_printf("The degree of %b is %b.\n",a[i],b[i]);
  }

  for (i=0;i<7;i++) {
    sollya_lib_clear_obj(a[i]);
    sollya_lib_clear_obj(b[i]);
  }

  sollya_lib_close();
  return 0;
}

