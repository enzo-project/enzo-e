#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

git diff -b | $dir/awk/diff-org.awk > diff.org

