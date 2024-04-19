#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Too few arguments"
    echo "Usage: ./test.sh <image-name>"
    exit 1
fi

# Change `python3` to `python` if your python3.x is alias to `python`
python3 utils.py img-to-bin "data/in-image/$1.jpg" "data/in-bin/$1.bin"

mpicc -o exec edge-detection.c && 
mpirun -np 8 exec "data/in-bin/$1.bin" "data/out-bin/$1.bin" &&
python3 utils.py bin-to-img "data/out-bin/$1.bin" "data/out-image/$1.jpg"