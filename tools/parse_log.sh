#!/usr/bin/awk -f

BEGIN{
    print "* TODO [/] Enzo-P / Cello hg log [C-u C-c # to update stats]"
}
{
    if (index($1,":") == length($1)) 
	{
#	    print NF;
	    if (NF > 0) {
		if ($1 == "changeset:") 
		    printf("** TODO %s", substr($1,1,length($1)-1));
		else if ($1 == "description:")
		    printf("*** %s\n", substr($1,1,length($1)-1));
		else if ($1 == "files:")
		    printf("*** TODO [/] %s", substr($1,1,length($1)-1));
		else
		    printf("*** %s",  substr($1,1,length($1)-1));
		}
	    if (NF > 1) {
		    printf ": ";
		    for (i=2; i<=NF; i++) {
			if ($1 == "files:") {
			    printf ("\n**** TODO [[file:%s][%s]]",$i,$i);
			    } else {
				printf ("%s ",$i); 
			    }
		    }
		    printf ("\n");
		}
	} else print
}
