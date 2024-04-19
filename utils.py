import sys
import logging
import argparse
from PIL import Image
from numpy import array, int32, uint8, fromfile

logging.basicConfig(level=logging.INFO)


def img_to_bin(img_filename, bin_filename):
    img = Image.open(img_filename).convert("L")
    width, height = img.size
    pixels = list(img.getdata())
    pixels = [pixels[i * width : (i + 1) * width] for i in range(height)]
    pixels = array(pixels, dtype=uint8)
    with open(bin_filename, "wb") as file:
        file.write(array([len(pixels[0]), len(pixels)], dtype=int32).tobytes())
        pixels.tofile(file)
    logging.info(f"Image {img_filename} converted to {bin_filename}")
    logging.info(f"Image size: {width} x {height}")


def bin_to_img(bin_filename, img_filename):
    pixels = array([])
    with open(bin_filename, "rb") as file:
        size = fromfile(file, dtype=int32, count=2)
        logging.debug(f"Image size: {size[0]} x {size[1]}")
        pixels = fromfile(file, dtype=uint8).reshape(size[::-1])
    img = Image.fromarray(pixels, "L")
    img.save(img_filename)
    logging.info(f"Binary {bin_filename} converted to {img_filename}")
    logging.info(f"Image size: {size[0]} x {size[1]}")


util_parser = argparse.ArgumentParser(description="Convert between image and text")
util_parser.add_argument(
    "mode", help="mode of operation", choices=["img-to-bin", "bin-to-img"]
)
util_parser.add_argument("filename", help="file name without path", nargs=2)

if __name__ == "__main__":
    args = util_parser.parse_args(sys.argv[1:])
    if args.mode == "img-to-bin":
        img_to_bin(args.filename[0], args.filename[1])
    elif args.mode == "bin-to-img":
        bin_to_img(args.filename[0], args.filename[1])
