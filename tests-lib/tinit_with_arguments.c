#include <sollya.h>

int callback(sollya_msg_t msg, void *data) {
  (void)msg; /* Avoiding "unused parameter" warning */
  (void)data; /* Avoiding "unused parameter" warning */
  return 0;
}

int main(void) {
  char *argv[] = { "5", "Hello world", "exp(x)", "", "Coucou" };
  int argc = 5;
  sollya_obj_t temp;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback, NULL);

  temp = sollya_lib_parse_string("__argv");
  sollya_lib_printf("__argv = %b\n", temp);
  sollya_lib_clear_obj(temp);

  sollya_lib_close();

  sollya_lib_init_with_arguments(argc, argv);
  sollya_lib_install_msg_callback(callback, NULL);

  temp = sollya_lib_parse_string("__argv");
  sollya_lib_printf("__argv = %b\n", temp);
  sollya_lib_clear_obj(temp);

  sollya_lib_close();
  return 0;
}
