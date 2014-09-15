#!/bin/bash

function mem {
   echo "$1 $2"
   if [ "$2" == "all" ]; then
     awk 'BEGIN{c=0}; /bytes-high/{print c,$6; c=c+1;}' < $1 > mem/mem-$2.data
   else
      awk 'BEGIN{c=0}; /Performance '"$2"' bytes-high/{print c,$6; c=c+1;}' < $1 > mem/mem-$2.data
   fi
}

if [ ! -d mem ]; then 
  mkdir mem
fi

mem $1 adapt
mem $1 refresh
mem $1 compute
mem $1 output
mem $1 stopping
mem $1 all

cd mem
$HOME/Cello/cello-src/util/ch-mem.py
