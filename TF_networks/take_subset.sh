#!/bin/bash

# use a very stupid method to take the subset of lines of $1 that contains the keywords in $2

while read line
do
    grep -w $line $1
done < $2
