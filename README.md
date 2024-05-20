# SF2568-project

## File structure

```
.
├── README.md               
├── data                # The data dir for script to use
│   ├── in-image        # The original images
│   ├── in-bin          # The bin files from original images
│   ├── out-image       # The bin files after edge detection
│   └── out-bin         # New images after edge detection
├── edge-detection.c    
├── utils.py            # Translating original images to bin file
├── test.sh             # Shell script for testing
└── requirements.txt    # pip install -r [this file]
```

## Pipeline Overview

Color Image -- (Python script) --> Gray Scale Imaging as 2D array in Binary File 


Binary File (original image) -- Edge detection (in C) --> Binary File (new image only contains edge)

Binary File  -- (Python script) --> Gray Scale Imaging (only contains edge)

## User guide
1. Prepare the Python environment, where Python is 3.10.x, and install the libraries using `pip install -r requirement.txt`. Prepare the C compiler, either mpicc or cc, configured with an MPI library.
2. First, prepare the image in any format as long as it is accepted by my Python library `Pillow`.
2. Use the following command to pre-process the image to the binary file.
    ```bash
    python3 utils.py img-to-bin [input image file] [output binary file]
    ```
3. Compile the source code `edge-detection.c` and put the option `-lm` at the end. For example, `cc edge-detection.c -O3 -lm`.
4. Run the executable. The `executable` takes two arguments: the input binary file and the output binary file. The example call on PDC will be,
    ```bash
    srun [executable] [input binary file] [output binary file]
    ```
5. Use the following command to post-process the binary file to the image.
    ```bash
    python3 utils.py img-to-bin [out binary file]  [input image file] [input image file]
    ```
    Note that the Python script always outputs a PNG image regardless of the argument's format to ensure lossless compression.

## Auto Threshold

By commenting out the line 
```C
#define AUTO_THRESHOLD 1
```
You can test the code without the auto threshold.

## Scripts
1. `test.sh`, used to speed up the test locally. To run this script, you MUST configure mpicc (from openmpi) and already put the corresponding image file in the directory `/data/in-image/`. This script used as follows

    ```bash
    ./test.sh [image file without path and postfix] [image file postfix]
    ```

    Example call: I have put Southbank-3.jpg in the directory as `/data/in-image/Southbank-3.jpg`. Using this script will give you the output image in `/data/out-image/Southbank-3.png`.

    ```bash
    ./test.sh Southbank-3 jpg
    ```