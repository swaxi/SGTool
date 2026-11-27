"""
Make a PowerPoint deck of every processed example in the SGTool reference.

Every example figure that docs/reference.md shows (docs/images/examples/*.png,
in the order of the manual, each once) is placed nine to a slide in a 3 x 3
grid, scaled to fit its cell. Each picture gets a small label, the figure's
alt text in the manual (e.g. "Low-latitude RTP"), at its top left, over the
input ("before") panel. The plain example grid (example_grid.png) is left out:
it is the input, not a processed example.

Needs python-pptx (pip install python-pptx). Run from anywhere:

    python scripts/make_examples_pptx.py                    # -> SGTool_examples.pptx
    python scripts/make_examples_pptx.py -o my_deck.pptx
"""

import argparse
import os
import re
import sys

try:
    from pptx import Presentation
    from pptx.dml.color import RGBColor
    from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
    from pptx.util import Emu, Inches, Pt
except ImportError:
    sys.exit("This script needs python-pptx: pip install python-pptx")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE = os.path.join(REPO, "docs", "reference.md")
SKIP = {"example_grid.png"}  # the input grid, not a processed example

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)  # 16:9
MARGIN, GAP = Inches(0.3), Inches(0.15)
COLS = ROWS = 3
LABEL_PT = 10


def examples():
    """(alt text, image path) of every example figure of the reference, in
    order, each once."""
    with open(REFERENCE, encoding="utf-8") as f:
        text = f.read()
    seen, out = set(), []
    for alt, rel in re.findall(r"!\[([^\]]*)\]\((images/examples/[^)\s]+)\)", text):
        name = os.path.basename(rel)
        if name in seen or name in SKIP:
            continue
        path = os.path.join(REPO, "docs", rel)
        if not os.path.exists(path):
            print(f"missing, skipped: {rel}")
            continue
        seen.add(name)
        out.append((alt.strip() or os.path.splitext(name)[0], path))
    return out


def image_size(path):
    """Pixel width and height of a PNG (read from its header, no Pillow needed)."""
    with open(path, "rb") as f:
        head = f.read(24)
    return int.from_bytes(head[16:20], "big"), int.from_bytes(head[20:24], "big")


def add_label(slide, text, left, top, max_width):
    """A small black-on-white text box at (left, top), as wide as its text."""
    width = min(max_width, Emu(int(Pt(LABEL_PT) * 0.6 * (len(text) + 1))))
    box = slide.shapes.add_textbox(left, top, width, Pt(LABEL_PT * 1.6))
    box.name = f"Label: {text}"
    box.fill.solid()
    box.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    tf = box.text_frame
    tf.margin_left = tf.margin_right = Pt(3)
    tf.margin_top = tf.margin_bottom = Pt(1)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.word_wrap = False
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = text
    run.font.size = Pt(LABEL_PT)
    run.font.bold = True
    run.font.name = "Arial"
    run.font.color.rgb = RGBColor(0, 0, 0)


def build(output):
    items = examples()
    prs = Presentation()
    prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
    blank = prs.slide_layouts[6]
    cell_w = (SLIDE_W - 2 * MARGIN - (COLS - 1) * GAP) // COLS
    cell_h = (SLIDE_H - 2 * MARGIN - (ROWS - 1) * GAP) // ROWS
    per_slide = COLS * ROWS
    for start in range(0, len(items), per_slide):
        slide = prs.slides.add_slide(blank)
        for k, (label, path) in enumerate(items[start:start + per_slide]):
            row, col = divmod(k, COLS)
            x0 = MARGIN + col * (cell_w + GAP)
            y0 = MARGIN + row * (cell_h + GAP)
            px_w, px_h = image_size(path)
            scale = min(cell_w / px_w, cell_h / px_h)  # fit the cell, keep the shape
            w, h = int(px_w * scale), int(px_h * scale)
            x = x0 + (cell_w - w) // 2
            y = y0 + (cell_h - h) // 2
            pic = slide.shapes.add_picture(path, x, y, w, h)
            pic.name = label
            # top left, over the input (left) panel; a long label may reach further
            add_label(slide, label, x, y, w)
    prs.save(output)
    print(f"{len(items)} examples on {len(prs.slides)} slides -> {output}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("-o", "--output", default="SGTool_examples.pptx",
                    help="the PowerPoint file to write (default: SGTool_examples.pptx)")
    build(ap.parse_args().output)
