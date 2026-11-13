#!/usr/bin/env python3
"""Compose four PyMOL renders into one publication figure."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analysis_results" / "generation_28_pymol"
OUTPUT = ROOT / "manuscript" / "figures" / "figure_generation_28_structures.png"
PANELS = [
    ("wild_type", "A)", "Wild type"),
    ("G978S", "B)", "G978S"),
    ("G984C", "C)", "G984C"),
    ("G987R", "D)", "G987R"),
]


def font(size: int, bold: bool = False):
    name = "Arial Bold.ttf" if bold else "Arial.ttf"
    path = Path("/System/Library/Fonts/Supplemental") / name
    return ImageFont.truetype(str(path), size=size)


def main() -> None:
    images = [Image.open(SOURCE / f"{name}.png").convert("RGB") for name, _, _ in PANELS]
    width = min(image.width for image in images)
    height = min(image.height for image in images)
    header = 92
    canvas = Image.new("RGB", (width * 2, (height + header) * 2), "white")
    draw = ImageDraw.Draw(canvas)
    for index, ((_, letter, title), image) in enumerate(zip(PANELS, images)):
        image = image.resize((width, height), Image.Resampling.LANCZOS)
        x = (index % 2) * width
        y = (index // 2) * (height + header)
        canvas.paste(image, (x, y + header))
        draw.text((x + 30, y + 18), letter, fill="black", font=font(46, bold=True))
        draw.text((x + 125, y + 20), title, fill="black", font=font(42, bold=True))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(OUTPUT, dpi=(300, 300), quality=95)


if __name__ == "__main__":
    main()
