#!/bin/bash

awk 'BEGIN {s=0}; {s=s+$2; print $1,s}' $1
