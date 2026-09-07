import argparse
from PIL import Image

p = argparse.ArgumentParser()
p.add_argument("input")
p.add_argument("output")
p.add_argument("--width", type=int, default=1000)
args = p.parse_args()

with Image.open(args.input) as image:
    ratio = args.width / image.width
    resized = image.convert("RGB").resize((args.width, round(image.height * ratio)), Image.Resampling.LANCZOS)
    resized.save(args.output, quality=82, optimize=True)
