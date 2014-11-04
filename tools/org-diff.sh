#!/bin/bash

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

hg diff | $dir/diff-org.awk > diff.org

