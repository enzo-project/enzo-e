#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: parse_error.sh < compiler-output > errors.org
#
# Generates org-mode file "errors.org" containing list of errors
# with links that can be opened with C-c C-e in emacs with org-mode
# installed.
#
# KNOWN ISSUES:
#
# * include/ files generating errors cannot be traced back to source
#   code directory, so it assumes src/Cello, which may be incorrect.
# * duplicate code with warning-org.awk
#----------------------------------------------------------------------


BEGIN {
    printf "* [[file:%s][ERRORS]] [/]\n",ARGV[1];
    c=0;
    ferror = 0;
}
/error:/ {
    c=c+1;
    split($1,a,":");
    b="";
    t=match ($1,"/");
    for (i=3; i<=NF; i++) {
        sub(/\[/,"(",$i);
        sub(/\]/,")",$i);
        b = (b " " $i);
    }
    print "** TODO [[file:"a[1]"::"a[2]"][Error "c": " b"]]"
}
/undefined reference to/ {
    c=c+1;
    split($1,a,":");
    b="";
    t=match ($1,"/");
    a[1] = "src"substr(a[1],t);
    sub(/include/,"src/Cello",a[1]);
    for (i=3; i<=NF; i++) {
        sub(/\[/,"(",$i);
        sub(/\]/,")",$i);
        b = (b " " $i);
    }
    print "** TODO [[file:"a[1]"::"a[2]"][Error "c": " b"]]"
}
{
    if (ferr == 1)
    {
        ferr = 0;
        c=c+1;
        split($1,a,":");
        b="";
        t=match ($1,"/");
        a[1] = "src"substr(a[1],t);
        sub(/include/,"src/Cello",a[1]);
        for (i=3; i<=NF; i++) {
            sub(/\[/,"(",$i);
            sub(/\]/,")",$i);
            b = (b " " $i);
        }
        split(a[2],A,".")
        print "** TODO [[file:"a[1]"::"A[1]"][Error "c": "ferror"]]"
    };
}
/Error:/ {
    ferr = 1;
    ferror = $0
    sub(/Error:/, "",ferror);
}
