#!/bin/bash

function mem {
   awk 'BEGIN{c=0}; /memory '"$2"'/{print c,$3; c=c+1;}' < $1 > mem/mem-$2.data
}

if [ ! -d mem ]; then 
  mkdir mem
fi

mem $1 adapt
mem $1 refresh
mem $1 compute
mem $1 output
mem $1 stopping

cd mem
$HOME/Cello/cello-src/util/ch-mem.py
