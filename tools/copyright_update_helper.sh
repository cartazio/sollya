#!/bin/sh

# Think also of updating sollya.tex, sollya.php (two places)
# and the meta copyright emitLegalNoticeAndDisclaimer in
# implement.c, as well as the corresponding reference files
# in the tests.

# This script explicitly exclude the following files, that must
# hence be manually checked and updated:
#    implement.c  because it contains emitLegalNoticeAndDisclaimer
#                 and there is therefore two notions of copyright
#                 for this file (the file itself, and the tool as
#                 a whole)
#    version.h    because the file itself is copyrighted, but it
#                 also contains the general copyright for the tool
#                 as a string.

# At the time of designing this script, only files of the form
# *.h *.c *.cpp *.y *.l are copyrighted. Also, no file in the
# subdirectories, except the above mentioned sollya.[tex|php] are
# copyrighted. But it might be good to check that it is still the
# case every time this script is run.

start_rev=sollya-4.1~
end_rev=HEAD

# A space separated list of hashes can be provided below. They will
# be excluded from the search. Typically useful when a particular
# commit changed a comma in all files of the repository
excluded_revs=""

# Methodology: the script outputs on stdout a list with one entry
# for each copyrighted file in Sollya. This entry first sums up the current
# copyright and contributors notice. Then, it lists all commits where the file
# has been modified since start_rev, together with the name of the contributor
# who performs the commit and the year of the commit.

# Note: we list years of the copyright as a range year1-year2, where year1 is
# the year of creation of the file, and year2 is the year of last modification,
# even if the file has not been modified during some year in-between.

cd /tmp
git clone https://scm.gforge.inria.fr/anonscm/git/sollya/sollya.git sollya-copyright
cd sollya-copyright

if [ -z "$excluded_revs" ]
then
    excluded_revs=ffffffffffffffffffffffffffffffffffffffff
else
    tmp_var=""
    for i in `printf "$excluded_revs" | sed -n 's/ /\n/g;p'`
    do
      tmp_var="$tmp_var"`git rev-parse $i`" "
    done
    excluded_revs=`printf "$tmp_var" | sed -n 's/ $//;p'`
fi

for file in *.c *.h *.cpp *.l *.y
do
  # WARNING: due to their particular content, these files
  # must be manually handled.
  if [ $file != "implement.c" -a $file != "version.h" ]
  then
    git log --full-history --pretty=format:'%H | %an | %ai' $start_rev..$end_rev $file | grep -v "^\("`printf "$excluded_revs" |sed -n 's/ /\\\\|/g;p'`"\)" | sed -n 's/^/   /;p' | sed -n 's/clauter/Christoph Lauter/;p' | sed -n 's/schevill/Sylvain Chevillard/;p' | sed -n 's/\(20[0-9][0-9]\)-.*$/\1/;p' > __sollya_tmp

    if grep "^[[:space:]]*[0-9a-f]" __sollya_tmp > /dev/null
    then
      notice=`grep -i Copyright $file`
      notice=$notice" "`grep -i Contributor $file`
      printf "%s" "$file:"
      printf "\n"
      printf "%s" "$notice"
      printf "\n"
      cat __sollya_tmp
      printf "\n\n"
    fi
    rm __sollya_tmp
  fi
done > copyright.log
mv copyright.log /tmp
cd /tmp
rm -rf sollya-copyright
cat copyright.log
rm copyright.log
