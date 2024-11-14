#!/bin/bash

if [ ! -f add_branches ]; then
    g++ -o add_branch add_branch.cpp `root-config --cflags --glibs`
fi

MAX_PROCESSES=10
running_processes=0

for sample in $(ls $1); do
    echo "Processing $sample"
    for f in $(ls $1/$sample/skimmed/*.root); do
        echo "Processing $f"
        ./add_branch $f &
        ((running_processes++))
        if ((running_processes >= MAX_PROCESSES)); then
            wait
            running_processes=0
        fi
    done
done

wait
