#include <sollya.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif

#if defined(HAVE_GETTIMEOFDAY) && HAVE_GETTIMEOFDAY

#include <inttypes.h>
#define LIST_SIZE 8000

int main(void) {
  sollya_obj_t listobj;
  sollya_obj_t tab[LIST_SIZE];
  sollya_obj_t a;
  int i, k;
  struct timeval t1;
  struct timeval t2;
  double length1;
  double length2;
  int64_t nb_iter;
  int res = 0;

  sollya_lib_init();

  for(i=0;i<LIST_SIZE;i++) tab[i] = SOLLYA_CONST(i);

  listobj = sollya_lib_list(tab, (LIST_SIZE/2));
  nb_iter = 1024;

  do {
    nb_iter = nb_iter*2;
    gettimeofday(&t1, NULL);
    for(k=0;k<nb_iter;k++) {
      sollya_lib_get_element_in_list(&a, listobj, (LIST_SIZE/2)-1);
      sollya_lib_clear_obj(a);
    }
    gettimeofday(&t2, NULL);
    length1 = (double)(t2.tv_sec - t1.tv_sec) * 1000000. +
      (double)(t2.tv_usec - t1.tv_usec);

  } while ( (nb_iter < INT64_MAX/2) && (length1 < 1000)); /* The loop must take at least 1000 ms */
  sollya_lib_clear_obj(listobj);

  /* Measures time on list 1 */
  listobj = sollya_lib_list(tab, (LIST_SIZE/2));
  gettimeofday(&t1, NULL);
  for(k=0;k<nb_iter;k++) {
    sollya_lib_get_element_in_list(&a, listobj, (LIST_SIZE/2)-1);
    sollya_lib_clear_obj(a);
  }
  gettimeofday(&t2, NULL);
  length1 = (double)(t2.tv_sec - t1.tv_sec) * 1000000. +
    (double)(t2.tv_usec - t1.tv_usec);
  sollya_lib_clear_obj(listobj);

  /* Measures time on list 2 */
  listobj = sollya_lib_list(tab, LIST_SIZE);

  listobj = sollya_lib_list(tab, 8000);
  gettimeofday(&t1, NULL);
  for(k=0;k<nb_iter;k++) {
    sollya_lib_get_element_in_list(&a, listobj, LIST_SIZE-1);
    sollya_lib_clear_obj(a);
  }
  gettimeofday(&t2, NULL);
  length2 = (double)(t2.tv_sec - t1.tv_sec) * 1000000. +
    (double)(t2.tv_usec - t1.tv_usec);
  sollya_lib_clear_obj(listobj);

  if ((0.75 <= length2/length1) && (length2/length1 <= 1.25))
    sollya_lib_printf("Testing that a call to get_element_in_list has complexity O(1): OK\n");
  else {
    sollya_lib_printf("Testing that a call to get_element_in_list has complexity O(1): not OK. Observed ratio: length1/length2 = %.0f/%.0f = %g (should be close to 1 for O(1) complexity, and close to 2 for linear complexity). Number of iterations: nb_iter = %" PRId64 "\n", length1, length2, length2/length1, nb_iter);
    res = 77;
  }
  
  for(i=0;i<LIST_SIZE;i++) sollya_lib_clear_obj(tab[i]);

  sollya_lib_close();

  return res;
}

#else /* gettimeofday not available: simply skip the test */
int main(void) { return 77; }
#endif
