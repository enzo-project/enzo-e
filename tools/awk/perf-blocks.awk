#!/usr/bin/awk -f

BEGIN{m=0};

/simulation num-blocks/ {
    l=substr($5,length($5));
    if (m<l) m=l;
    if (l==0) {
	for (i=0; i<=m; i++) {
	    printf "%g ",b[i];
	}
	printf "\n";
    }
    b[l]=$6
}
