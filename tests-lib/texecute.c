#include <sollya.h>

int callback(sollya_msg_t msg, void *data) {
  int message = sollya_lib_get_msg_id(msg);
  (void)data; /* Avoiding "unused parameter" warning */

  sollya_lib_printf("Unexpected message\n");
  return 0;
}

int main(void) {

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback, NULL);

  /* TODO */

  sollya_lib_close();
  return 0;
}
