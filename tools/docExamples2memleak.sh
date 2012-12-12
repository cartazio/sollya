#!/bin/sh

if [ $# -eq 0 ]
  then dir="."
  else dir=$1
fi

dir_bak=`pwd`
cd $dir

if ! /bin/ls *.shlp 1> /dev/null 2>&1
then
  echo "Usage $0 dir where dir is the directory containing shlp files"
  exit 1
fi

touch ___toto && rm ___toto
for i in *.shlp
do
  grep '#EXAMPLE' -A1000 $i >> ___toto
done
grep -v "SEEALSO" ___toto | tr -d '\n' | sed -n 's/#EXAMPLE/\n/g;p' > ___titi
mv ___titi ___toto

i=1
while [ $i -le `wc -l ___toto |sed -n 's/[^0-9]//g;p'` ]
do
  head -n $i ___toto |tail -n 1 > $dir_bak/$i.memleak
  i=`expr $i + 1`
done
rm ___toto

cd $dir_bak