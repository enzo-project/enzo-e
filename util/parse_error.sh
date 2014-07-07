#!/bin/bash

# usage: parse_error.sh <compiler-output>
#
# output: errors.org with links to code

echo "* [[file:test/out.scons][ERRORS]] [/]" > errors.org
awk 'BEGIN{c=0}; /error:/{c=c+1; split($1,a,":"); b="";  sub(/build/, "src",a[1]); sub(/include/,"src/Cello",a[1]); for (i=3; i<=NF; i++) { sub(/\[/,"(",$i); sub(/\]/,")",$i); b = (b " " $i); } print "** TODO [[file:"a[1]"::"a[2]"][Error "c": " b"]]"}' test/out.scons >> errors.org

echo "* [[file:test/out.scons][WARNINGS]] [/]" >> errors.org
awk 'BEGIN{c=0}; /warning:/{c=c+1; split($1,a,":"); b="";  sub(/build/, "src",a[1]); sub(/include/,"src/Cello",a[1]); for (i=3; i<=NF; i++) { sub(/\[/,"(",$i); sub(/\]/,")",$i); b = (b " " $i); } print "** TODO [[file:"a[1]"::"a[2]"][Warning "c": " b"]]"}' test/out.scons >> errors.org


