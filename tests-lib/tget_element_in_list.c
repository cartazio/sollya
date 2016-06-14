#include <sollya.h>
#include <inttypes.h>

typedef struct __sollya_time_struct_t sollya_time_t;
int sollya_gettime(sollya_time_t *);
int64_t sollya_timediff_ms(sollya_time_t *, sollya_time_t *);
sollya_time_t *sollya_gettime_var();

#define LIST_SIZE 8000

int main(void) {
  sollya_obj_t listobj;
  sollya_obj_t tab[LIST_SIZE];
  sollya_obj_t a;
  int i, k;
  sollya_time_t *t1;
  sollya_time_t *t2;
  int64_t length1;
  int64_t length2;
  int64_t nb_iter;

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


  /* Time Complexity */
  for(i=0;i<LIST_SIZE;i++) tab[i] = SOLLYA_CONST(i);

  t1 = sollya_gettime_var();
  t2 = sollya_gettime_var();

  listobj = sollya_lib_list(tab, (LIST_SIZE/2));
  nb_iter = 1024;
  do {
    nb_iter = nb_iter*2;
    sollya_gettime(t1);
    for(k=0;k<nb_iter;k++) {
      sollya_lib_get_element_in_list(&a, listobj, (LIST_SIZE/2)-1);
      sollya_lib_clear_obj(a);
    }
    sollya_gettime(t2);
    length1 = sollya_timediff_ms(t1,t2);
  } while ( (nb_iter < INT64_MAX/2) && (length1 < 1000)); /* The loop must take at least 1000 ms */
  sollya_lib_clear_obj(listobj);

  /* Measures time on list 1 */
  listobj = sollya_lib_list(tab, (LIST_SIZE/2));
  sollya_gettime(t1);
  for(k=0;k<nb_iter;k++) {
    sollya_lib_get_element_in_list(&a, listobj, (LIST_SIZE/2)-1);
    sollya_lib_clear_obj(a);
  }
  sollya_gettime(t2);
  length1 = sollya_timediff_ms(t1,t2);
  sollya_lib_clear_obj(listobj);

  /* Measures time on list 2 */
  listobj = sollya_lib_list(tab, LIST_SIZE);
  sollya_gettime(t1);
  for(k=0;k<nb_iter;k++) {
    sollya_lib_get_element_in_list(&a, listobj, LIST_SIZE-1);
    sollya_lib_clear_obj(a);
  }
  sollya_gettime(t2);
  length2 = sollya_timediff_ms(t1,t2);
  sollya_lib_clear_obj(listobj);
  sollya_lib_free(t1);
  sollya_lib_free(t2);

  if ((0.75 <= ((double)length2)/((double)length1)) && (((double)length2)/((double)length1) <= 1.25))
    sollya_lib_printf("Testing that a call to get_element_in_list has complexity O(1): OK\n");
  else
    sollya_lib_printf("Testing that a call to get_element_in_list has complexity O(1): not OK. Observed ratio: length1/length2 = %" PRId64 "/%" PRId64 " = %g (should be close to 1 for O(1) complexity, and close to 2 for linear complexity). Number of iterations: nb_iter = %" PRId64 "\n", length1, length2, ((double)length2)/((double)length1), nb_iter);

  for(i=0;i<LIST_SIZE;i++) sollya_lib_clear_obj(tab[i]);

 sollya_lib_close();

  return 0;
}
