#include <sollya.h>

int main(void) {
  sollya_obj_t listobj;
  sollya_obj_t a;
  int i;

  sollya_lib_init();

  /* Empty list */
  listobj = sollya_lib_build_list(NULL);
  a = NULL;
  i = 0;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* Negative index on a regular list */
  listobj = sollya_lib_build_list(SOLLYA_CONST(3.),
                                  SOLLYA_COS(SOLLYA_X_),
                                  sollya_lib_double_obj(),
                                  NULL);
  a = NULL;
  i = -1;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* Out-of-bounds index on a regular list */
  listobj = sollya_lib_build_list(SOLLYA_CONST(3.),
                                  SOLLYA_COS(SOLLYA_X_),
                                  sollya_lib_double_obj(),
                                  NULL);
  a = NULL;
  i = 3;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* First element on a regular list */
  listobj = sollya_lib_build_list(SOLLYA_CONST(3.),
                                  SOLLYA_COS(SOLLYA_X_),
                                  sollya_lib_double_obj(),
                                  NULL);
  a = NULL;
  i=0;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* Last element on a regular list */
  listobj = sollya_lib_build_list(SOLLYA_CONST(3.),
                                  SOLLYA_COS(SOLLYA_X_),
                                  sollya_lib_double_obj(),
                                  NULL);
  a = NULL;
  i=2;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* Negative index on an end-elliptic list */
  listobj = sollya_lib_build_end_elliptic_list(SOLLYA_CONST(3.),
                                               SOLLYA_COS(SOLLYA_X_),
                                               sollya_lib_double_obj(),
                                               NULL);
  a = NULL;
  i=-1;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* index after the ellipsis on an end-elliptic list ending with a non-integer element */
  listobj = sollya_lib_build_end_elliptic_list(SOLLYA_CONST(3.),
                                               SOLLYA_COS(SOLLYA_X_),
                                               sollya_lib_double_obj(),
                                               NULL);
  a = NULL;
  i=3;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* index after the ellipsis on an end-elliptic list ending with an integer element */
  listobj = sollya_lib_build_end_elliptic_list(SOLLYA_PI,
                                               SOLLYA_COS(SOLLYA_X_),
                                               SOLLYA_CONST(4.),
                                               NULL);
  a = NULL;
  i=3;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* First element on an end-elliptic list */
  listobj = sollya_lib_build_end_elliptic_list(SOLLYA_PI,
                                               SOLLYA_COS(SOLLYA_X_),
                                               SOLLYA_CONST(4.),
                                               NULL);
  a = NULL;
  i=0;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* Last element before the ellipsis on an end-elliptic list */
  listobj = sollya_lib_build_end_elliptic_list(SOLLYA_PI,
                                               SOLLYA_COS(SOLLYA_X_),
                                               SOLLYA_CONST(4.),
                                               NULL);
  a = NULL;
  i=2;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);


  /* Not a list */
  listobj = SOLLYA_COS(SOLLYA_X_);
  a = NULL;
  i=0;
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Success on %b : %b (i=%d)\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n", listobj, i, a);
  }
  if (sollya_lib_get_element_in_list(&a, listobj, i)) {
    sollya_lib_printf("Second call on %b : %b (i=%d)\n\n", listobj, a, i);
    sollya_lib_clear_obj(a);
  }
  else {
    sollya_lib_printf("Failure on %b (i=%d and a=%b).\n\n", listobj, i, a);
  }
  sollya_lib_clear_obj(listobj);

 sollya_lib_close();

  return 0;
}
