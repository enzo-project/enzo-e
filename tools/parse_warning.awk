#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: parse_warning.sh < compiler-output > warnings.org
#
# Generates org-mode file "warnings.org" containing list of warnings
# with links that can be opened with C-c C-e in emacs with org-mode
# installed.
#
# KNOWN ISSUES:
#
# * include/ files generating warnings cannot be traced back to source
#   code directory, so it assumes src/Cello, which may be incorrect.
# * duplicate code with parse_error.sh
#----------------------------------------------------------------------


BEGIN {
    print "* [[file:test/out.scons][WARNINGS]] [/]";
    c=0;
    fwarning = 0;
}
/warning:/ {
    c=c+1; 
    split($1,a,":"); 
    b="";  
    sub(/build/, "src",a[1]); 
    sub(/include/,"src/Cello",a[1]);
    for (i=3; i<=NF; i++) {
	sub(/\[/,"(",$i);
	sub(/\]/,")",$i);
	b = (b " " $i); 
    }
    print "** TODO [[file:"a[1]"::"a[2]"][Warning "c": " b"]]"
}
{
    if (ferr == 1) 
	{
	    ferr = 0;
	    c=c+1; 
	    split($1,a,":"); 
	    b="";  
	    sub(/build/, "src",a[1]); 
	    sub(/include/,"src/Cello",a[1]);
	    for (i=3; i<=NF; i++) {
		sub(/\[/,"(",$i);
		sub(/\]/,")",$i);
		b = (b " " $i); 
	    }
	    split(a[2],A,".")
	    print "** TODO [[file:"a[1]"::"A[1]"][Warning "c": "fwarning"]]"
	};
	
}
/Warning:/ {
    ferr = 1;
    fwarning = $0
    sub(/Warning:/, "",fwarning); 
}
