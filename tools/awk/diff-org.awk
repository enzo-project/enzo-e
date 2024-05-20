#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: org-mode-diff.sh < <"git diff" output > diff.org
#
#----------------------------------------------------------------------


BEGIN{
    p=1;
    l=1
}
/^diff --git/ { # git
   file=$NF; p=0;
   sub("b/","",file);
   n=split(file,s,"/");
   shortfile = s[n];
}
/^\-\-\-/ {p=0}
/^index /{p=0}
/^\+\+\+/ {p=0}
/^@@/ {
    line=$3;
    line=substr(line,2,index(line,",")-2);
    print "*** TODO [[file:" file"::"line"]["l " " shortfile ":" line "]]"
    lp = l
};
{
    if (p==1)
	print;
    p=1
}
{
    l=l+1
}
END {print l}


