#!/bin/bash

~/.cellorc/msg-counts.awk < $1 > msg-counts.data
gnuplot ~/.cellorc/msg-counts.gnu
