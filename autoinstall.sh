#!/bin/sh

GMP_PKG_NAME=gmp-5.1.2
MPFR_PKG_NAME=mpfr-3.1.2
MPFI_PKG_NAME=mpfi-1.5.1
FPLLL_PKG_NAME=libfplll-4.0.3
SOLLYA_PKG_NAME=sollya-4.0

LOG_DIR=$HOME/log
PREFIX_PATH=$HOME/local
SOLLYA_CFLAGS="-g -g3 -gdwarf-2 -O0 -Wall -Wextra -ansi"
SOLLYA_CXXFLAGS="-g -g3 -gdwarf-2 -O0 -Wall -Wextra -ansi"

###############################################################################


if which wget 1> /dev/null 2>&1
then
  GET_PRG="wget --no-check-certificate"
else
  if which curl 1> /dev/null 2>&1
  then
    GET_PRG="curl -O"
  else
    printf "Please provide either wget or curl to use this script\n"
    exit 1
  fi
fi

printf "\nDownloading up-to-date autoinstall script... "
$GET_PRG http://www-sop.inria.fr/members/Sylvain.Chevillard/autoinstall-core 1> /dev/null 2>&1
if test "$?" -ne 0
then
  printf "Failed.\n"
  exit 1
fi
printf "done.\n"

. ./autoinstall-core
rm autoinstall-core