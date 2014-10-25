#!/usr/bin/awk -f

#----------------------------------------------------------------------
#
# usage: ls -1 <file-glob> | parse_ls.sh > ls.org
#
# Generates org-mode file "ls.org" containing TODO links to files.
#
#----------------------------------------------------------------------


BEGIN{
    print "* TODO Files [/] [C-u C-c # to update stats]"
}
{
    file = $1;
    print "** TODO [[file:" file"::"file"]["file"]]";
}
