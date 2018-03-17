#!/usr/bin/awk -f


BEGIN{
    p=0;
    print "* TODO [/] [C-u C-c # to update stats]"
}
/lost in loss record/ {
    p=1;
    printf "** TODO %d/%d %s\n",$(NF-2),$NF,$2
}

/\(/ && /\)/ && /pp:/ {
    if (p==1) {
	n=length($NF);
	file=substr($NF,2,n-2);
	sub(/:/,"::",file);
	printf "   - [ ] [[file:source/%s][%s]]\n",file,file
    }
}
