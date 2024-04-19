#!/bin/bash

# Change `python3` to `python` if your python3.x is alias to `python`
# echo "Type in the filename of the image, not path and postfix are needed"
python3 utils.py img-to-bin "data/in-image/$1.jpg" "data/in-bin/$1.bin"

mpicc -o exec edge-detection.c && 
mpirun -np 8 exec "data/in-bin/$1.bin" "data/out-bin/$1.bin" &&
## p data/in-bin/$1.bin data/out-bin/$1.bin &&
python3 utils.py bin-to-img "data/out-bin/$1.bin" "data/out-image/$1.jpg"