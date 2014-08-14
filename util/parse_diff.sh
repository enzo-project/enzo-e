#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: parse_diff.sh < <"hg diff" output > diff.org
#
# Generates org-mode file "diff.org" containing list of warnings
# with links that can be opened with C-c C-e in emacs with org-mode
# installed.
#
# KNOWN ISSUES:
#
# * include/ files generating warnings cannot be traced back to source
#   code directory, so it assumes src/Cello, which may be incorrect.
# * duplicate code with parse_error.sh
#----------------------------------------------------------------------


BEGIN{
    p=1;
    print "* TODO [/] [C-u C-c # to update stats]"
}
/diff -r/ {
   file=$4; p=0;
   n=split(file,s,"/");
   shortfile = s[n];
   print "** TODO [/]",shortfile;
}
/^\-\-\-/ {p=0}
/^\+\+\+/ {p=0}
/@@/ {
    line=$3;
    line=substr(line,2,index(line,",")-2);
    print "*** TODO [[file:" file"::"line"]["shortfile ":" line "]]"};
{
    if (p==1)
	print;
    p=1
}
