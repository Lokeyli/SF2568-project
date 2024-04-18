# SF2568-project

## File structure

```
.
├── README.md               
├── data
│   ├── in-image        # The original images
│   ├── in-txt          # The text files from original images
│   ├── out-image       # The text files after edge detection
│   └── out-txt         # New images after edge detection
├── edge-detection.c    
├── imgToTxt.py         # Translating original images to text file
└── requirements.txt    # pip install -r [this file]
```

## Pipeline Overview

Color Image -- (Python script) --> Gray Scale Imaging as 2D array in Text File 


Text File (original image) -- Edge detection (in C) --> Text File (new image only contains edge)

Text File  -- (Python script) --> Gray Scale Imaging (only contains edge)

## Pipeline Detail
The translation from color image to gray-gradient image represented in as txt file is done by python script `imgToTxt.py`, by calling below in the command line

`python imgToTxt.py [image file path] [txt file path]`