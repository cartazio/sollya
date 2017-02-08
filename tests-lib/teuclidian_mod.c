#include <sollya.h>

static const char* inputs[] = { "1001", 
				"231", 
				"13", 
				"17", 
				"-210", 
				"462", 
				"6/7", 
				"33/13", 
				"exp(13)", 
				"sin(17)", 
				"24 + 68 * _x_ + 74 * _x_^2 + 39 * _x_^3 + 10 * _x_^4 + _x_^5", 
				"480 + 776 * _x_ + 476 * _x_^2 + 138 * _x_^3 + 19 * _x_^4 + _x_^5", 
				"1001 * _x_^2", 
				"231 * _x_", 
				"exp(_x_)", 
				"_x_^2", 
				"-14", 
				"15", 
				"-213", 
				"-5", 
				"23/13", 
				"11/17", 
				"exp(13)", 
				"-sin(17)", 
				"4 + 4 * _x_ + _x_^2", 
				"2 * _x_^3", 
				"_x_^3", 
				"sin(_x_)", 
				NULL };

int main(void) {
  int i, j, nb;
  sollya_obj_t *a;
  sollya_obj_t r;

  sollya_lib_init();

  for (nb=0; inputs[nb] != NULL; nb++);
  
  a = (sollya_obj_t *) sollya_lib_calloc(nb, sizeof(sollya_obj_t));
  for (i=0; i<nb; i++) a[i] = sollya_lib_parse_string(inputs[i]);

  for (i=0; i<nb; i++) {
    for (j=0; j<nb; j++) {
      r = sollya_lib_euclidian_mod(a[i], a[j]);
      sollya_lib_printf("mod(%b, %b) = %b\n", a[i], a[j], r);
      sollya_lib_clear_obj(r);
    }
  }
  
  for (i=0; i<nb; i++) sollya_lib_clear_obj(a[i]);
  sollya_lib_free(a);
  
  sollya_lib_close();
  return 0;
}

