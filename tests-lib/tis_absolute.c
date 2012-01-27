#include <sollya.h>

int main(void) {
  sollya_obj_t a[51];
  int t[51];
  int i;

  sollya_lib_init();

  a[0] = sollya_lib_absolute();
  a[1] = sollya_lib_binary();
  a[2] = sollya_lib_decimal();
  a[3] = sollya_lib_default();
  a[4] = sollya_lib_dyadic();
  a[5] = sollya_lib_error();
  a[6] = sollya_lib_false();
  a[7] = sollya_lib_file();
  a[8] = sollya_lib_fixed();
  a[9] = sollya_lib_floating();
  a[10] = sollya_lib_hexadecimal();
  a[11] = sollya_lib_honorcoeffprec();
  a[12] = sollya_lib_off();
  a[13] = sollya_lib_on();
  a[14] = sollya_lib_perturb();
  a[15] = sollya_lib_pi();
  a[16] = sollya_lib_postscript();
  a[17] = sollya_lib_postscriptfile();
  a[18] = sollya_lib_powers();
  a[19] = sollya_lib_relative();
  a[20] = sollya_lib_true();
  a[21] = sollya_lib_void();
  a[22] = sollya_lib_constant_from_double(1.0);
  a[23] = sollya_lib_free_variable();
  a[24] = sollya_lib_string("absolute");
  a[25] = sollya_lib_string("");
  a[26] = sollya_lib_list(NULL, 0);
  a[27] = sollya_lib_list(a, 2);
  a[28] = sollya_lib_end_elliptic_list(a, 2);
  a[29] = sollya_lib_range(sollya_lib_constant_from_int(1),
                           sollya_lib_constant_from_int(2));
  a[30] = sollya_lib_round_down();
  a[31] = sollya_lib_round_up();
  a[32] = sollya_lib_round_towards_zero();
  a[33] = sollya_lib_round_to_nearest();
  a[34] = sollya_lib_double_obj();
  a[35] = sollya_lib_build_function_double(sollya_lib_free_variable());
  a[36] = sollya_lib_single_obj();
  a[37] = sollya_lib_build_function_single(sollya_lib_free_variable());
  a[38] = sollya_lib_quad_obj();
  a[39] = sollya_lib_build_function_quad(sollya_lib_free_variable());
  a[40] = sollya_lib_halfprecision_obj();
  a[41] = sollya_lib_build_function_halfprecision(sollya_lib_free_variable());
  a[42] = sollya_lib_doubleextended_obj();
  a[43] = sollya_lib_build_function_doubleextended(sollya_lib_free_variable());
  a[44] = sollya_lib_double_double_obj();
  a[45] = sollya_lib_build_function_double_double(sollya_lib_free_variable());
  a[46] = sollya_lib_triple_double_obj();
  a[47] = sollya_lib_build_function_triple_double(sollya_lib_free_variable());
  a[48] = sollya_lib_build_function_cos(sollya_lib_constant_from_int(5));
  a[49] = sollya_lib_end_elliptic_list(NULL, 0);
  a[50] = sollya_lib_parse_string("{.a = 5; .b = exp(_x_)};");

  for(i=0;i<51;i++)  t[i] = sollya_lib_is_absolute(a[i]);

  if (t[0]) sollya_lib_printf("sollya_lib_is_absolute detects absolute\n");
  if (t[1]) sollya_lib_printf("sollya_lib_is_absolute detects binary\n");
  if (t[2]) sollya_lib_printf("sollya_lib_is_absolute detects decimal\n");
  if (t[3]) sollya_lib_printf("sollya_lib_is_absolute detects default\n");
  if (t[4]) sollya_lib_printf("sollya_lib_is_absolute detects dyadic\n");
  if (t[5]) sollya_lib_printf("sollya_lib_is_absolute detects error\n");
  if (t[6]) sollya_lib_printf("sollya_lib_is_absolute detects false\n");
  if (t[7]) sollya_lib_printf("sollya_lib_is_absolute detects file\n");
  if (t[8]) sollya_lib_printf("sollya_lib_is_absolute detects fixed\n");
  if (t[9]) sollya_lib_printf("sollya_lib_is_absolute detects floating\n");
  if (t[10]) sollya_lib_printf("sollya_lib_is_absolute detects hexadecimal\n");
  if (t[11]) sollya_lib_printf("sollya_lib_is_absolute detects honorcoeffprec\n");
  if (t[12]) sollya_lib_printf("sollya_lib_is_absolute detects off\n");
  if (t[13]) sollya_lib_printf("sollya_lib_is_absolute detects on\n");
  if (t[14]) sollya_lib_printf("sollya_lib_is_absolute detects perturb\n");
  if (t[15]) sollya_lib_printf("sollya_lib_is_absolute detects pi\n");
  if (t[16]) sollya_lib_printf("sollya_lib_is_absolute detects postscript\n");
  if (t[17]) sollya_lib_printf("sollya_lib_is_absolute detects postscriptfile\n");
  if (t[18]) sollya_lib_printf("sollya_lib_is_absolute detects powers\n");
  if (t[19]) sollya_lib_printf("sollya_lib_is_absolute detects relative\n");
  if (t[20]) sollya_lib_printf("sollya_lib_is_absolute detects true\n");
  if (t[21]) sollya_lib_printf("sollya_lib_is_absolute detects void\n");
  if (t[22]) sollya_lib_printf("sollya_lib_is_absolute detects a constant\n");
  if (t[23]) sollya_lib_printf("sollya_lib_is_absolute detects free_variable\n");
  if (t[24]) sollya_lib_printf("sollya_lib_is_absolute detects the string \"absolute\"\n");
  if (t[25]) sollya_lib_printf("sollya_lib_is_absolute detects the empty string.\n");
  if (t[26]) sollya_lib_printf("sollya_lib_is_absolute detects the empty list.\n");
  if (t[27]) sollya_lib_printf("sollya_lib_is_absolute detects a non empty list.\n");
  if (t[28]) sollya_lib_printf("sollya_lib_is_absolute detects a non empty end-elliptic list.\n");
  if (t[29]) sollya_lib_printf("sollya_lib_is_absolute detects a range.\n");
  if (t[30]) sollya_lib_printf("sollya_lib_is_absolute detects RD.\n");
  if (t[31]) sollya_lib_printf("sollya_lib_is_absolute detects RU.\n");
  if (t[32]) sollya_lib_printf("sollya_lib_is_absolute detects RZ.\n");
  if (t[33]) sollya_lib_printf("sollya_lib_is_absolute detects RN.\n");
  if (t[34]) sollya_lib_printf("sollya_lib_is_absolute detects D as symbol.\n");
  if (t[35]) sollya_lib_printf("sollya_lib_is_absolute detects D(_x_).\n");
  if (t[36]) sollya_lib_printf("sollya_lib_is_absolute detects S as symbol.\n");
  if (t[37]) sollya_lib_printf("sollya_lib_is_absolute detects S(_x_).\n");
  if (t[38]) sollya_lib_printf("sollya_lib_is_absolute detects QD as symbol.\n");
  if (t[39]) sollya_lib_printf("sollya_lib_is_absolute detects QD(_x_).\n");
  if (t[40]) sollya_lib_printf("sollya_lib_is_absolute detects HP as symbol.\n");
  if (t[41]) sollya_lib_printf("sollya_lib_is_absolute detects HP(_x_).\n");
  if (t[42]) sollya_lib_printf("sollya_lib_is_absolute detects DE as symbol.\n");
  if (t[43]) sollya_lib_printf("sollya_lib_is_absolute detects DE(_x_).\n");
  if (t[44]) sollya_lib_printf("sollya_lib_is_absolute detects DD as symbol.\n");
  if (t[45]) sollya_lib_printf("sollya_lib_is_absolute detects DD(_x_).\n");
  if (t[46]) sollya_lib_printf("sollya_lib_is_absolute detects TD as symbol.\n");
  if (t[47]) sollya_lib_printf("sollya_lib_is_absolute detects TD(_x_).\n");
  if (t[48]) sollya_lib_printf("sollya_lib_is_absolute detects a constant expression.\n");
  if (t[49]) sollya_lib_printf("sollya_lib_is_absolute detects an empty end-elliptic list.\n");
  if (t[50]) sollya_lib_printf("sollya_lib_obj_is_structure detects structure { .a = 5; .b = exp(_x_)}.\n");

  for(i=0;i<51;i++) sollya_lib_clear_obj(a[i]);
  sollya_lib_close();
  return 0;
}
