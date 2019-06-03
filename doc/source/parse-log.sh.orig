#!/bin/csh -f

awk 'BEGIN {p=0}; /changeset:   /{p=0; print "*revision",substr($2,1,4),"*"}; {if (p==1) print} /description:/{p=1}; ' $1
