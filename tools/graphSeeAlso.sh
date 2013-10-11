#!/bin/sh

# Small script to be run in the doc directory: it scans all the .shlp files,
# finds the #SEEALSO directives and constructs a file graph.dot describing the
# graph. Then, it renders the graph with the dot tool.

printf "Generating the graph...\n"
printf "digraph seealsos {\n" > graph.dot

grep "=" keywords.def | \
sed -n 's/^\([^=]*\).*/\1/;p' | \
while read ligne
do
 printf "$ligne [label=\"$ligne\"]\n"
done >> graph.dot


for i in *.shlp
do
 j=`printf $i | sed -n 's/\.shlp//;p' | tr 'a-z' 'A-Z'`
 grep SEEALSO $i | sed -n 's/#SEEALSO $\(.*\)/'$j'->\1/;p'
done >> graph.dot

printf "}" >> graph.dot

printf "Rendering the graph...\n"
dot graph.dot -Tpng -o graph.png