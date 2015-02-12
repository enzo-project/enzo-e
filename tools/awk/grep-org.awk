#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: grep -n <text> <file-glob> | grep-org.awk > grep.org
#
# Generates org-mode file "grep.org" containing TODO links to grep
# output. 
#
#----------------------------------------------------------------------


BEGIN{
    print "* TODO [/] [C-u C-c # to update stats]"
}
{
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
