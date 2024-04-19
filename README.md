# SF2568-project

## File structure

```
.
├── README.md               
├── data
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

## Pipeline Detail
The translation between color image and gray-gradient image represented in as binary file is done by python script `utils.py`. Use the following command to see the detail

```bash
python utils.py -h
```

## Scripts
1. `test.sh`, used as followed

    ```bash
    ./test.sh [image file without path and postfix]
    ```

    This is intended to speed up the test, assuming you have already put the corresponding `.JPEG` file in the directory `/data/in-image/`.