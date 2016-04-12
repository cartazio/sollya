#include "sollya.h"

int cmp(int *res, void **args) {
  sollya_obj_t *true_args = (sollya_obj_t *)args;
  *res = sollya_lib_cmp_objs_structurally(true_args[0], true_args[1]);
  return 1;
}

int hcmp(int *res, void **args) {
  sollya_obj_t *true_args = (sollya_obj_t *)args;
  int64_t h1,h2;
  h1 = sollya_lib_hash(true_args[0]);
  h2 = sollya_lib_hash(true_args[1]);
  *res = (h1 == h2);
  return 1;
}

