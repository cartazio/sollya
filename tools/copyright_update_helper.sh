############# IMPORTANT NOTE ###############
# This script is not up to date anymore since when we switched from svn to git to manage Sollya
# It must be updated.

#!/bin/sh

start_rev=2023
end_rev=HEAD
excluded_revs=""

cd /tmp
svn co svn+ssh://schevill@scm.gforge.inria.fr/svn/sollya/trunk sollya-copyright
cd sollya-copyright

if [ -z "$excluded_revs" ]
  then excluded_revs=`expr $start_rev - 1`
fi

for file in `grep --recursive -i -l "Copyright.*20" * | grep -v "\.svn"`
do
  svn log $file -r$start_rev:$end_rev | grep "^r[0-9]" | sed -n 's/^\(r[0-9]* | [^ ]* | [0-9]*\).*/\1/;p' | grep -v `printf "$excluded_revs" |sed -n 's/ /\\\\|r/g;p' |sed -n 's/^/r/;p'` | sed -n 's/^/   /;p' > __sollya_tmp
  if grep "^[[:space:]]*r[0-9]" __sollya_tmp > /dev/null
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
done > copyright.log

# Be aware that there could be false positives. Sometimes, a commit does not require to change the copyright notice (typically, a change of syntax over the whole source code). Such revisions can be explicitely excluded by adding them to excluded_revs variable.
# Also, be aware that the work done in branches is not reflected by copyright.log.

# Methodology: the file /tmp/sollya-copyright/copyright.log contains an entry for each copyrighted file in Sollya. This entry first sums up the current copyright and contributors notice. Then, it lists all commits where the file has been modified since start_rev, together with the name of the contributor who performs the commit and the year of the commit.

# Note: we list years of the copyright as a range year1-year2, where year1 is the year of creation of the file, and year2 is the year of last modification, even if the file has not been modified during some year in-between.

