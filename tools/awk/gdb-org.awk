#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: gdb-org.sh < gdb frame > gdb.org
#
# Generates org-mode file "gdb.org" containing stack frame with links
# that can be opened with C-c C-o in emacs with org-mode installed.
#
#----------------------------------------------------------------------


BEGIN {
    print "* GDB stack frame [/]";
    c=0;
    ferror = 0;
}
/^#/ {
    enum=$1
}
/ at / {
    n=length($0)
    i=index($0," at ")
    s = substr($0,i+4,n-i)
    t=s
    i=index(t,"/")
    while (i != 0) {
	# remove all chars left of right-most '/'
	n=length(t)
	t=substr(t,i+1,n-i);
	i=index(t,"/")
    }
    gsub(/:/,"::",s)
    gsub(/build/,"src",s)
    
    print "** TODO",enum,"[[file:" s "][" t "]]"
}
