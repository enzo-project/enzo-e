#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: doxygen ... | doxygen-org.awk > doxygen.org
#
# Generates org-mode file "doxygen.org" containing TODO links to doxygen
# warnings. 
#
#----------------------------------------------------------------------


BEGIN{
    print "* TODO [/] [C-u C-c # to update stats]"
}
/warning: /{
    line=$0;
    file_end=index(line,":")
    file = substr(line,1,file_end-1)
    line = substr(line,file_end+1,length(line)-file_end+1)
    num_end=index(line,":");
    num  = substr(line,1,num_end-1)
    line = substr(line,num_end+1,length(line)-num_end+1)
#    print file,num,line
    print "** TODO [[file:" file"::"num"]["file":"num"]]"

};
