#!/usr/bin/awk -f

BEGIN{
    p=0;
    print "* TODO [/] [C-u C-c # to update stats]"
}
/lost in loss record/ {
    p=1;
    for (i=1; i<NF; i++) {
	if ($i == "blocks") nb=$(i-1);
    }
    printf "** TODO %d/%d %s bytes in %s blocks\n",$(NF-2),$NF,$2,nb
}
/Conditional jump or move/ {
    p=1;
    for (i=1; i<NF; i++) {
	if ($i == "blocks") nb=$(i-1);
    }
    printf "** TODO Invalid %s of size %d\n",$3,$NF
}    
/Invalid/ {
    p=1;
    for (i=1; i<NF; i++) {
	if ($i == "blocks") nb=$(i-1);
    }
    printf "** TODO Conditional jump or move depends on uninitialised value(s)\n"
}
/Use of uninitialised value of size/ {
    p=1;
    for (i=1; i<NF; i++) {
	if ($i == "blocks") nb=$(i-1);
    }
    printf "** TODO Use of uninitialised value of size %d\n",$8
}
/\(/ && /\)/ && /pp:/ || /F:/ {
    if (p==1) {
	n=length($NF);
	file=substr($NF,2,n-2);
	sub(/:/,"::",file);
	printf "   - [ ] [[file:src/%s][%s]]\n",file,file
    }
}
