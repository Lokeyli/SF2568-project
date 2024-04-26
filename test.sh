#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Too few arguments, should be 2"
    echo "Usage: ./test.sh <image-name> <image-format>"
    exit 1
fi

# Change `python3` to `python` if your python3.x is alias to `python`
python3 utils.py img-to-bin "data/in-image/$1.$2" "data/in-bin/$1.bin"

mpicc -o exec edge-detection.c -lm && 
mpirun -np 8 exec "data/in-bin/$1.bin" "data/out-bin/$1.bin" &&
python3 utils.py bin-to-img "data/out-bin/$1.bin" "data/out-image/$1.$2"