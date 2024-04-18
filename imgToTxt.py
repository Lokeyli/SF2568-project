import sys
from PIL import Image


def img_to_txt(img_filename, txt_filename):
    img = Image.open(img_filename).convert("L")
    width, height = img.size
    pixels = list(img.getdata())
    pixels = [pixels[i * width : (i + 1) * width] for i in range(height)]
    with open(txt_filename, "w") as file:
        for row in pixels:
            file.write(" ".join(str(p) for p in row) + "\n")


if __name__ == "__main__":
    args = sys.argv[1:]
    img_to_txt(args[0], args[1])
