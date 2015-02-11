#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: org-mode-diff.sh < <"hg diff" output > diff.org
#
#----------------------------------------------------------------------


BEGIN{
    p=1;
    print "* TODO [/] [C-u C-c # to update stats]"
}
/diff -r/ {
   file=$NF; p=0;
   n=split(file,s,"/");
   shortfile = s[n];
}
/^\-\-\-/ {p=0}
/^\+\+\+/ {p=0}
/@@/ {
    line=$3;
    line=substr(line,2,index(line,",")-2);
    print "** TODO [[file:" file"::"line"]["shortfile ":" line "]]"};
{
    if (p==1)
	print;
    p=1
}
