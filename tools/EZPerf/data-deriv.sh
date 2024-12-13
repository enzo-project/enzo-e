#!/bin/bash

awk 'BEGIN {p=0}; {print $1,$2-p; p=$2}' $1
