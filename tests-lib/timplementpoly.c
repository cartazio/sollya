#include <sollya.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>


#define A_DIM 15
#define TEMP_IMPLEMENTATION_FILE "timplementpoly.implementation.c"
#define TEMP_PROOF_FILE "timplementpoly.implementation.gappa"
#define READBUFFERSIZE 1024

int callback(int message) {
  switch(message) {
  case SOLLYA_MSG_IMPLEMENTED_POLY_DIFFERS_FROM_ORIGINAL_ONE:
    sollya_lib_printf("Caught the message: the implemented polynomial is different from the original one.\n");
    break;
  case SOLLYA_MSG_INFERED_COEFF_PREC_HIGHER_THAN_REQUIRED:
    sollya_lib_printf("Caught the message: the infered precision of a coefficient is higher than what seems to be needed to meet the accuracy target.\n");
    break;    
  case SOLLYA_MSG_COEFF_NOT_TWICE_GREATER_THAN_SUBPOLY:
    sollya_lib_printf("Caught the message: a coefficient in a Horner scheme is not guaranteed to also be twice as large as the subpolynomial.\n");
    break;    
  case SOLLYA_MSG_ERROR_ON_DETERMINING_THE_REQUIRED_PRECISIONS:
    sollya_lib_printf("Caught the message: an error has occured during the determination of the required precisions.\n");
    break;    
  default:
    sollya_lib_printf("Unexpected warning %d.\n", message);
  }
  return 0;
}

void read_and_print_file(char *filename) {
  FILE *fd;
  char readBuffer[READBUFFERSIZE];
  int i;
  size_t readChars;

  if ((fd = fopen(filename,"r")) == NULL) {
    exit(1);
  }

  while (1) {
    for (i=0;i<READBUFFERSIZE;i++) {
      readBuffer[i] = '\0';
    }
    readChars = fread(readBuffer, sizeof(char), READBUFFERSIZE - 1, fd);
    for (i=0;i<readChars;i++) {
      if (readBuffer[i] == '\0') readBuffer[i] = ' ';
    }
    printf("%s",readBuffer);
    if (readChars < (READBUFFERSIZE - 1) * sizeof(char)) break;
  }
  printf("\n");

  if (fclose(fd)) {
    exit(1);
  }

}

int main(void) {
  sollya_obj_t a[A_DIM];
  sollya_obj_t res, temp;
  int i;
  struct stat statBuffer;
  int status;
  int errornum;
  int toRemove;

  sollya_lib_init();
  sollya_lib_install_msg_callback(callback);

  sollya_lib_set_display(temp = sollya_lib_dyadic());
  sollya_lib_clear_obj(temp);


  /* implementpoly(1 - TD(1/6) * x^2,[-1b-10;1b-10],1b-60,DD,"p","implementation.c"); */
  for (i=0;i<A_DIM;i++) {
    a[i] = NULL;
  }

  toRemove = 1;
  status = stat(TEMP_IMPLEMENTATION_FILE, &statBuffer);
  if (status != 0) {
    errornum = errno;
    if (errornum != ENOENT) {
      sollya_lib_close();
      return 1;
    }
    toRemove = 0;
  }
  if (toRemove) {
    if (remove(TEMP_IMPLEMENTATION_FILE) != 0) {
      sollya_lib_close();
      return 1;
    }
  }

  a[0] = sollya_lib_parse_string("1 - TD(1/6) * x^2");
  a[1] = sollya_lib_parse_string("[-1b-10;1b-10]");
  a[2] = sollya_lib_parse_string("1b-60");
  a[3] = sollya_lib_double_double_obj();
  a[4] = sollya_lib_string("p");
  a[5] = sollya_lib_string(TEMP_IMPLEMENTATION_FILE);
  
  res = sollya_lib_implementpoly(a[0],a[1],a[2],a[3],a[4],a[5],NULL);

  sollya_lib_printf("implementpoly(%b,%b,%b,%b,%b,\"%b\") returns %b and produces the following code:\n\n",
		    a[0],a[1],a[2],a[3],a[4],a[5],res);

  read_and_print_file(TEMP_IMPLEMENTATION_FILE);
  
  for (i=0;i<A_DIM;i++) {
    if (a[i] != NULL) {
      sollya_lib_clear_obj(a[i]);
    }
  }
  sollya_lib_clear_obj(res);

  /* implementpoly(1 - simplify(TD(1/6)) * x^2,[-1b-10;1b-10],1b-60,DD,"p","implementation.c",honorcoeffprec); */
  for (i=0;i<A_DIM;i++) {
    a[i] = NULL;
  }

  toRemove = 1;
  status = stat(TEMP_IMPLEMENTATION_FILE, &statBuffer);
  if (status != 0) {
    errornum = errno;
    if (errornum != ENOENT) {
      sollya_lib_close();
      return 1;
    }
    toRemove = 0;
  }
  if (toRemove) {
    if (remove(TEMP_IMPLEMENTATION_FILE) != 0) {
      sollya_lib_close();
      return 1;
    }
  }

  a[0] = sollya_lib_parse_string("1 - TD(1/6) * x^2");
  a[1] = sollya_lib_parse_string("[-1b-10;1b-10]");
  a[2] = sollya_lib_parse_string("1b-60");
  a[3] = sollya_lib_double_double_obj();
  a[4] = sollya_lib_string("p");
  a[5] = sollya_lib_string(TEMP_IMPLEMENTATION_FILE);
  a[6] = sollya_lib_honorcoeffprec();
  
  res = sollya_lib_implementpoly(a[0],a[1],a[2],a[3],a[4],a[5],a[6],NULL);

  sollya_lib_printf("implementpoly(%b,%b,%b,%b,%b,\"%b\",%b) returns %b and produces the following code:\n\n",
		    a[0],a[1],a[2],a[3],a[4],a[5],a[6],res);

  read_and_print_file(TEMP_IMPLEMENTATION_FILE);
  
  for (i=0;i<A_DIM;i++) {
    if (a[i] != NULL) {
      sollya_lib_clear_obj(a[i]);
    }
  }
  sollya_lib_clear_obj(res);

  /* implementpoly(0x3ff0000000000000 + x * (0x3ff0000000000000 + x * (0x3fe0000000000000 + x * (0x3fc5555555555559 + x * (0x3fa55555555555bd + x * (0x3f811111111106e2 + x * (0x3f56c16c16bf5eb7 + x * (0x3f2a01a01a292dcd + x * (0x3efa01a0218a016a + x * (0x3ec71de360331aad + x * (0x3e927e42e3823bf3 + x * (0x3e5ae6b2710c2c9a + x * (0x3e2203730c0a7c1d + x * 0x3de5da557e0781df)))))))))))),[-1/2;1/2],1b-60,D,"p","implementation.c",honorcoeffprec,"implementation.gappa"); */
  for (i=0;i<A_DIM;i++) {
    a[i] = NULL;
  }

  toRemove = 1;
  status = stat(TEMP_IMPLEMENTATION_FILE, &statBuffer);
  if (status != 0) {
    errornum = errno;
    if (errornum != ENOENT) {
      sollya_lib_close();
      return 1;
    }
    toRemove = 0;
  }
  if (toRemove) {
    if (remove(TEMP_IMPLEMENTATION_FILE) != 0) {
      sollya_lib_close();
      return 1;
    }
  }

  toRemove = 1;
  status = stat(TEMP_PROOF_FILE, &statBuffer);
  if (status != 0) {
    errornum = errno;
    if (errornum != ENOENT) {
      sollya_lib_close();
      return 1;
    }
    toRemove = 0;
  }
  if (toRemove) {
    if (remove(TEMP_PROOF_FILE) != 0) {
      sollya_lib_close();
      return 1;
    }
  }

  a[0] = sollya_lib_parse_string("0x3ff0000000000000 + x * (0x3ff0000000000000 + x * (0x3fe0000000000000 + x * (0x3fc5555555555559 + x * (0x3fa55555555555bd + x * (0x3f811111111106e2 + x * (0x3f56c16c16bf5eb7 + x * (0x3f2a01a01a292dcd + x * (0x3efa01a0218a016a + x * (0x3ec71de360331aad + x * (0x3e927e42e3823bf3 + x * (0x3e5ae6b2710c2c9a + x * (0x3e2203730c0a7c1d + x * 0x3de5da557e0781df))))))))))))");
  a[1] = sollya_lib_parse_string("[-1/2;1/2]");
  a[2] = sollya_lib_parse_string("1b-60");
  a[3] = sollya_lib_double_double_obj();
  a[4] = sollya_lib_string("p");
  a[5] = sollya_lib_string(TEMP_IMPLEMENTATION_FILE);
  a[6] = sollya_lib_honorcoeffprec();
  a[7] = sollya_lib_string(TEMP_PROOF_FILE);
  
  res = sollya_lib_implementpoly(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],NULL);

  sollya_lib_printf("implementpoly(%b,%b,%b,%b,%b,\"%b\",%b,\"%b\") returns %b and produces the following code:\n\n",
		    a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res);

  read_and_print_file(TEMP_IMPLEMENTATION_FILE);

  sollya_lib_printf("Additionnally, the following Gappa proof is produced:\n\n");

  read_and_print_file(TEMP_PROOF_FILE);
  
  for (i=0;i<A_DIM;i++) {
    if (a[i] != NULL) {
      sollya_lib_clear_obj(a[i]);
    }
  }
  sollya_lib_clear_obj(res);

  toRemove = 1;
  status = stat(TEMP_IMPLEMENTATION_FILE, &statBuffer);
  if (status != 0) {
    errornum = errno;
    if (errornum != ENOENT) {
      sollya_lib_close();
      return 1;
    }
    toRemove = 0;
  }
  if (toRemove) {
    if (remove(TEMP_IMPLEMENTATION_FILE) != 0) {
      sollya_lib_close();
      return 1;
    }
  }

  toRemove = 1;
  status = stat(TEMP_PROOF_FILE, &statBuffer);
  if (status != 0) {
    errornum = errno;
    if (errornum != ENOENT) {
      sollya_lib_close();
      return 1;
    }
    toRemove = 0;
  }
  if (toRemove) {
    if (remove(TEMP_PROOF_FILE) != 0) {
      sollya_lib_close();
      return 1;
    }
  }

  sollya_lib_close();
  return 0;
}
