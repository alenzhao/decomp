#!/bin/bash

# use a very stupid method to take the subset of lines of $1 that contains the keywords in $2

while read line
do
    # take the first matching line
    grep -m 1 -w -P "\t${line}\t" $1
done < $2
