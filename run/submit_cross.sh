#!/bin/bash
for i in 0 1 2 3 4
do
    echo "Cross Val fold $i"
    ./cross_val.sh $i &
done

