#include <sollya.h>
#include <mpfr.h>
#include <mpfi.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(void) {
  sollya_obj_t f[13];
  int i, j;
  uint64_t h1, h2;

  sollya_lib_init();

  f[0] = SOLLYA_CONST_SI64(17);
  f[1] = SOLLYA_CONST(17.0);
  f[2] = SOLLYA_CONST_SI64(42);
  f[3] = SOLLYA_CONST_SI64(25);
  f[4] = sollya_lib_sub(f[2],f[3]);
  f[5] = sollya_lib_parse_string("1 + 17 * x + 136 * x^2 + 680 * x^3 + 2380 * x^4 + 6188 * x^5 + 12376 * x^6 + 19448 * x^7 + 24310 * x^8 + 24310 * x^9 + 19448 * x^10 + 12376 * x^11 + 6188 * x^12 + 2380 * x^13 + 680 * x^14 + 136 * x^15 + 17 * x^16 + x^17");
  f[6] = sollya_lib_parse_string("(1 + x)^17");
  f[7] = sollya_lib_parse_string("(1 + x)^42");
  f[8] = sollya_lib_degree(f[6]);
  f[9] = sollya_lib_parse_string("proc (a,b) { \"When I add \", b, \" to \", a, \", I get \", (a + b); return (a + b); }");
  f[10] = SOLLYA_CONST_SI64(8);
  f[11] = SOLLYA_CONST_SI64(9);
  f[12] = sollya_lib_apply(f[9], f[10], f[11], NULL);

  sollya_lib_printf("Objects:\n");
  for (i=0;i<=12;i++) {
    sollya_lib_printf("%02d: %b\n", i, f[i]);
  }
  sollya_lib_printf("\n");
  
  sollya_lib_printf("Hash equality table:\n");
  sollya_lib_printf("-- ");
  for (j=0;j<=12;j++) {
    sollya_lib_printf("%02d ",j);
  }
  sollya_lib_printf("\n");
  for (i=0;i<=12;i++) {
    h1 = sollya_lib_hash(f[i]);
    sollya_lib_printf("%02d ", i);
    for (j=0;j<=12;j++) {
      h2 = sollya_lib_hash(f[j]);
      if (h1 == h2) {
	sollya_lib_printf(" * ");
      } else {
	sollya_lib_printf("   ");
      }
    }
    sollya_lib_printf("\n");
  }

  for(i=0;i<=12;i++) sollya_lib_clear_obj(f[i]);
  
  sollya_lib_close();
  
  return 0;
}
