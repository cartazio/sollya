#include <sollya.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

#define FILENAME_SUFFIX ".tmp"

int callback(sollya_msg_t msg, void *data) {
  (void)msg; /* Avoiding "unused parameter" warning */
  (void)data; /* Avoiding "unused parameter" warning */
  return 0;
}

int main(int argc, char **argv) {
  char *tempFile;
  int fopenerr, fcloseerr;
  FILE *fd;
  sollya_obj_t fileObj;

  if (argc < 1) return 1;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback, NULL);

  tempFile = (char *) sollya_lib_calloc(strlen(argv[0]) + strlen(FILENAME_SUFFIX) + 1, sizeof(char));
  strcpy(tempFile, argv[0]);
  strcpy(tempFile + strlen(argv[0]), FILENAME_SUFFIX);
  if ((fd = fopen(tempFile, "w")) == NULL) {
    fopenerr = errno;
    fprintf(stderr, "Error opening temporary file \"%s\": \"%s\"\n", tempFile, strerror(fopenerr));
    sollya_lib_free(tempFile);
    sollya_lib_close();
    return 1;
  }
  fprintf(fd, "f = exp(x);\nf(5);\nprocedure a(f,n) {\n    \"f = \", f;\n    \"n = \", n;\n    f(n);\n    return f(n);\n};\nt = a(f + sin(x),5);\nprint(t);\n");
  if (fclose(fd) != 0) {
    fcloseerr = errno;
    fprintf(stderr, "Error closing temporary file \"%s\": \"%s\"\n", tempFile, strerror(fcloseerr));
    sollya_lib_free(tempFile);
    sollya_lib_close();
    return 1;
  }
  fileObj = sollya_lib_string(tempFile);
  sollya_lib_free(tempFile);
  sollya_lib_execute(fileObj);
  sollya_lib_clear_obj(fileObj);
 
  sollya_lib_close();
  return 0;
}

