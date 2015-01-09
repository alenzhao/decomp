#!/bin/bash

# the list of file names is in $@

cat <(cut -f1 "$@") <(cut -f2 "$@") | sort -u
