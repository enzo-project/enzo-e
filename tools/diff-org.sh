#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

hg diff -b | $dir/awk/diff-org.awk > diff.org

