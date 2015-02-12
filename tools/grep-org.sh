#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

grep -n $1 $2 | $dir/awk/grep-org.awk > grep.org

