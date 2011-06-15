#/bin/bash

beforedir=`pwd`
cd $1
echo "#include \"expansion.h\"" > impl.c
echo " " >> impl.c
cat $beforedir/$2  >> impl.c
gcc -Winline -std=c99 -D$3 -DPOLYNOMIALNAME=$4 -fPIC  -c impl.c 
gcc -Winline -std=c99 -D$3 -DPOLYNOMIALNAME=$4 -fPIC  -c expansion.c 
gcc -std=c99 -shared -o $5 impl.o expansion.o -lgmp -lmpfr 
#rm impl.c
rm impl.o
rm expansion.o
cd $beforedir


