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
/Invalid/ {
    p=1;
    for (i=1; i<NF; i++) {
	if ($i == "blocks") nb=$(i-1);
    }
    printf "** TODO Invalid %s of size %d\n",$(NF-4),$NF
}

/\(/ && /\)/ && /pp:/ || /F:/ {
    if (p==1) {
	n=length($NF);
	file=substr($NF,2,n-2);
	sub(/:/,"::",file);
	printf "   - [ ] [[file:source/%s][%s]]\n",file,file
    }
}
