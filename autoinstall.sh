#!/bin/sh

# Usage: if you wish to use a particular version of GMP, MPFR, MPFI,
# FPLLL or Sollya, set the following variables to the name of the
# package. If the variable is left empty, the standard version
# available on the system will be used.
#
# The special package name 'sollya-git' is available for downloading sollya
# from the git repository.
#
# Since fplll is not distributed as a package anymore, the special package
# names 'libfplll-git' and 'libfplll-git-somehash' are also available to use
# respectively the git master or some particular git commit (identified by its
# hash 'somehash') of fplll's Github
GMP_PKG_NAME=gmp-5.1.3
MPFR_PKG_NAME=mpfr-3.1.2
MPFI_PKG_NAME=mpfi-1.5.1
FPLLL_PKG_NAME=libfplll-4.0.4
SOLLYA_PKG_NAME=sollya-4.1

# If $LOG_DIR is empty, the current directory will be used.  If
# $PREFIX_PATH is empty, everything will be installed in the standard
# place of the system
LOG_DIR=$HOME/log
PREFIX_PATH=$HOME/local

GMP_CONFIGURE_OPTIONS=""
MPFR_CONFIGURE_OPTIONS=""
MPFI_CONFIGURE_OPTIONS=""
FPLLL_CONFIGURE_OPTIONS=""

SOLLYA_CFLAGS="-g -g3 -gdwarf-2 -O0 -Wall -Wextra -std=c99"
SOLLYA_CXXFLAGS="-g -g3 -gdwarf-2 -O0 -Wall -Wextra"

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
rm -rf autoinstall-core

$GET_PRG http://www-sop.inria.fr/members/Sylvain.Chevillard/autoinstall-core 1> /dev/null 2>&1
if test "$?" -ne 0
then
  printf "Failed.\n"
  exit 2
fi
printf "done.\n"

. ./autoinstall-core
rm -rf autoinstall-core
