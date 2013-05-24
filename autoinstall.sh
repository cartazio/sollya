#!/bin/sh

GMP_PKG_NAME=gmp-5.1.2
MPFR_PKG_NAME=mpfr-3.1.2
MPFI_PKG_NAME=mpfi-1.5.1
FPLLL_PKG_NAME=libfplll-4.0.3
SOLLYA_PKG_NAME=sollya-svn

LOG_DIR=$HOME/log
PREFIX_PATH=$HOME/.local
SOLLYA_CFLAGS="-g -g3 -gdwarf-2 -O0 -Wall -Wextra -std=c99 -ansi"
SOLLYA_CXXFLAGS="-g -g3 -gdwarf-2 -O0 -Wall -Wextra -ansi"

###############################################################################

failure_check() {
 if test "$?" -ne 0
 then
   printf "Failed.\n"
   exit 1
 fi
}

if test -z "$LOG_DIR"
then
  LOG_DIR=$HOME/log
fi
if test ! -d "$LOG_DIR"
then
  printf "Creating $LOG_DIR... "
  mkdir "$LOG_DIR"; failure_check
  printf "done.\n"
fi

if test -n "$PREFIX_PATH" -a ! -d "$PREFIX_PATH"
then
  printf "Creating $PREFIX_PATH... "
  mkdir "$PREFIX_PATH"; failure_check
  printf "done.\n"
fi

if test -n "$PREFIX_PATH"
then
  PREFIX="--prefix=$PREFIX_PATH"
  WITHGMP="--with-gmp=$PREFIX_PATH"
  WITHMPFR="--with-mpfr=$PREFIX_PATH"
  WITHMPFI="--with-mpfi=$PREFIX_PATH"
  WITHFPLLL="--with-fplll=$PREFIX_PATH"
else
  PREFIX=""
  WITHGMP=""
  WITHMPFR=""
  WITHMPFI=""
  WITHFPLLL=""
fi

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

######################### INSTALLATION OF GMP #########################
if test -n "$GMP_PKG_NAME"
then
  printf "\nDownloading GMP... "
  $GET_PRG ftp://ftp.gmplib.org/pub/$GMP_PKG_NAME/$GMP_PKG_NAME.tar.bz2 1> /dev/null 2>&1 ; failure_check
  printf "done.\n"

  tar -jxvf $GMP_PKG_NAME.tar.bz2 > /dev/null
  cd $GMP_PKG_NAME

  printf "Running configure... "
  ./configure "$PREFIX" > "$LOG_DIR/gmp_configure.log" ; failure_check
  printf "done.\n"

  printf "Compiling... "
  make > "$LOG_DIR/gmp_compile.log" ; failure_check
  printf "done.\n"

  printf "Running make check... "
  make check > "$LOG_DIR/gmp_check.log" ; failure_check
  printf "done.\n"

  printf "Installing it... "
  make install > "$LOG_DIR/gmp_install.log" ; failure_check
  printf "done.\n"

  cd ..
  rm -rf $GMP_PKG_NAME.tar.bz2
  rm -rf $GMP_PKG_NAME
  printf "\n"
fi

######################### INSTALLATION OF MPFR ########################
if test -n "$MPFR_PKG_NAME"
then
  printf "Downloading MPFR... "
  $GET_PRG http://www.mpfr.org/mpfr-current/$MPFR_PKG_NAME.tar.bz2 1> /dev/null 2>&1 ; failure_check
  printf "done.\n"

  tar -jxvf $MPFR_PKG_NAME.tar.bz2 > /dev/null
  cd $MPFR_PKG_NAME

  printf "Running configure... "
  ./configure $PREFIX $WITHGMP > "$LOG_DIR/mpfr_configure.log" ; failure_check
  printf "done.\n"

  printf "Compiling... "
  make > "$LOG_DIR/mpfr_compile.log" ; failure_check
  printf "done.\n"

  printf "Running make check... "
  make check > "$LOG_DIR/mpfr_check.log" ; failure_check
  printf "done.\n"

  printf "Installing it... "
  make install > "$LOG_DIR/mpfr_install.log" ; failure_check
  printf "done.\n"

  cd ..
  rm -rf $MPFR_PKG_NAME.tar.bz2
  rm -rf $MPFR_PKG_NAME
  printf "\n"
fi

######################### INSTALLATION OF MPFI ########################
if test -n "$MPFI_PKG_NAME"
then
  printf "Downloading MPFI... "
  $GET_PRG https://gforge.inria.fr/frs/download.php/30129/$MPFI_PKG_NAME.tar.bz2 1> /dev/null 2>&1 ; failure_check
  printf "done.\n"

  tar -jxvf $MPFI_PKG_NAME.tar.bz2 > /dev/null
  cd $MPFI_PKG_NAME

  printf "Running configure... "
  ./configure $PREFIX $WITHGMP $WITHMPFR > "$LOG_DIR/mpfi_configure.log" ; failure_check
  printf "done.\n"

  printf "Compiling... "
  make > "$LOG_DIR/mpfi_compile.log" ; failure_check
  printf "done.\n"

  printf "Running make check... "
  make check > "$LOG_DIR/mpfi_check.log" ; failure_check
  printf "done.\n"

  printf "Installing it... "
  make install > "$LOG_DIR/mpfi_install.log" ; failure_check
  printf "done.\n"

  cd ..
  rm -rf $MPFI_PKG_NAME.tar.bz2
  rm -rf $MPFI_PKG_NAME
  printf "\n"
fi

######################## INSTALLATION OF FPLLL ########################
if test -n "$FPLLL_PKG_NAME"
then
  printf "Downloading FPLLL... "
  $GET_PRG http://perso.ens-lyon.fr/damien.stehle/fplll/$FPLLL_PKG_NAME.tar.gz 1> /dev/null 2>&1 ; failure_check
  printf "done.\n"

  tar -zxvf $FPLLL_PKG_NAME.tar.gz > /dev/null
  cd $FPLLL_PKG_NAME

  if test -e fplll.patch
  then
    printf "Patching bug in FPLLL configure file... "
    patch -p0 < fplll.patch 1> /dev/null 2>&1 ; failure_check
    autoreconf -i 1> /dev/null 2>&1 ; failure_check
    printf "done.\n"
  fi
  printf "Running configure... "
  ./configure $PREFIX $WITHGMP $WITHMPFR > "$LOG_DIR/fplll_configure.log" ; failure_check
  printf "done.\n"

  printf "Compiling... "
  make > "$LOG_DIR/fplll_compile.log" ; failure_check
  printf "done.\n"

  printf "Running make check... "
  make check > "$LOG_DIR/fplll_check.log" ; failure_check
  printf "done.\n"

  printf "Installing it... "
  make install > "$LOG_DIR/fplll_install.log" ; failure_check
  printf "done.\n"

  cd ..
  rm -rf $FPLLL_PKG_NAME.tar.gz
  rm -rf $FPLLL_PKG_NAME
fi

######################## INSTALLATION OF SOLLYA ########################
printf "Downloading Sollya... "
case "$SOLLYA_PKG_NAME" in
  sollya-svn )
    svn co svn://scm.gforge.inria.fr/svn/sollya/trunk $SOLLYA_PKG_NAME > /dev/null ; failure_check
    ;;
  sollya-3.0 )
    SOLLYA_URL=https://gforge.inria.fr/frs/download.php/28570/sollya-3.0.tar.bz2
    ;;
  sollya-2.9 )
    SOLLYA_URL=https://gforge.inria.fr/frs/download.php/28249/sollya-2.9.tar.bz2
    ;;
  sollya-2.0 )
    SOLLYA_URL=https://gforge.inria.fr/frs/download.php/26857/sollya-2.0.tar.bz2
    ;;
  sollya-weekly-* )
    SOLLYA_URL=http://sollya.gforge.inria.fr/$SOLLYA_PKG_NAME.tar.bz2
    ;;
  * )
    printf "Failed. (Unknown package)\n"; exit 1
    ;;
  esac
if test "$SOLLYA_PKG_NAME" != "sollya-svn"
then
  $GET_PRG $SOLLYA_URL 1> /dev/null 2>&1 ; failure_check
fi
printf "done.\n"

if test "$SOLLYA_PKG_NAME" != "sollya-svn"
then
  tar -jxvf $SOLLYA_PKG_NAME.tar.bz2 > /dev/null
  cd "$SOLLYA_PKG_NAME"
else
  cd "$SOLLYA_PKG_NAME"
  printf "Generating configure file... "
  ./autogen.sh 1> /dev/null 2>&1 ; failure_check
  printf "done.\n"
fi

printf "Running configure... "
./configure $PREFIX $WITHGMP $WITHMPFR $WITHMPFI $WITHFPLLL CFLAGS="$SOLLYA_CFLAGS" CXX_FLAGS="$SOLLYA_CXXFLAGS" > "$LOG_DIR/sollya_configure.log" ; failure_check
printf "done.\n"

printf "Compiling... "
make > "$LOG_DIR/sollya_compile.log" ; failure_check
printf "done.\n"

printf "Running make check... "
make check > "$LOG_DIR/sollya_check.log" ; failure_check
printf "done.\n"

printf "Installing it... "
make install > "$LOG_DIR/sollya_install.log" ; failure_check
printf "done.\n"

cd ..


