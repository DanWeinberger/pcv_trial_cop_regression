"""
Extract data from the Andrejko 2024 Figure 2 (Vaccine 42:3559): "Cases included
in analysis stratified by serotype and receipt of PCV13".

Unlike the Whitney 2006 figure (which is vector geometry, see extract_figure_data.py),
Figure 2 here is embedded as a RASTER image, so the pdf2plot vector method does not
apply. Instead we extract the PNG and do colour-based pixel analysis:

  * The three stacked segments use ggplot2's default 3-colour palette:
        blue  (97,156,255)  = "3+ PCV13 doses"     (bottom of stack)
        green (0,186,56)     = "1-2 PCV13 doses"    (middle)
        red   (248,118,109)  = "No PCV doses"       (top)
  * Two panels: A) Early implementation (baseline row 278, y=80 at row 42)
                B) Late  implementation (baseline row 693, y=80 at row 457)
    Both scale at 236 px / 80 cases = 2.95 px per case (from major gridlines).
  * Bars sit on a regular x-grid: 44 serotypes, first (19A) center x=94,
    last (NT) center x=1858.
  * For each bar we take a central 30-px window and, per colour, divide the
    colored-pixel count by the number of bar columns to get the mean segment
    height, then convert to cases.

Output is in long/tidy format: one row per (period, serotype, vaccine status),
with Total_Cases giving the study-wide denominator (from Table 1) for that
vaccine status, for computing per-serotype rates.

Validation against the paper's text (all near-exact):
  * serotype 19A total (both panels) = 96  -> "43.0% (96/223) of VT-IPD"
  * serotype 19F: Early ~4 / Late ~41      -> "3.8% [early]... 34.4% [late]"
  * 3 + 19F totals ~= 105                   -> 19A+3+19F = 201/223
  * panel sums ~528 / ~647 vs Fig 1's 524 / 637; combined ~1170 vs 1161

Requires: pymupdf, pillow, numpy   (pip install pymupdf pillow numpy)
"""

import csv
from pathlib import Path

import fitz  # PyMuPDF
import numpy as np
from PIL import Image

# ggplot2 default 3-colour palette for the stacked segments
BLUE = (97, 156, 255)    # 3+ PCV13 doses
GREEN = (0, 186, 56)     # 1-2 PCV13 doses
RED = (248, 118, 109)    # No PCV doses
COLOR_TOL = 60

# y-axis calibration (image pixel rows), per panel
PANEL_A = dict(y_scan=(30, 279), baseline=278, top=42)   # Early
PANEL_B = dict(y_scan=(445, 694), baseline=693, top=457)  # Late
PX_PER_UNIT = 236 / 80.0   # (baseline - top) / 80 cases = 2.95

# x-grid of bar centres
X_FIRST, X_LAST = 94.0, 1858.0

# serotypes in plotting order (VT/PCV13-type first, then NVT)
SEROTYPES = [
    "19A", "3", "7F", "19F", "14", "9V", "18C", "23F", "6B", "6C",   # PCV13-type
    "22F", "33F", "15C", "35B", "15B", "38", "23B", "10A", "12F", "23A",
    "15A", "21", "11A", "35F", "8", "9N", "16F", "15B/C", "7C", "31",
    "20", "17F", "28A", "34", "11B", "10F", "33A", "15F", "18F", "35D",
    "24F/24A/24", "24F/A/B", "35B:35D", "NT",
]
# Study-wide denominators from Table 1 (VT cases + NVT controls), not derived
# from the figure. No doses of any PCV = 47 + 53 = 100. >=3 doses PCV13 =
# 108 + 600 = 708. 1-2 doses = (>=1 dose total 176+885=1061) - 708 = 353,
# which matches the text ("excluded 353 participants... 1 or 2 doses").
STUDY_TOTAL_NO_PCV = 100
STUDY_TOTAL_1_2_DOSES = 353
STUDY_TOTAL_3PLUS_DOSES = 708


def color_mask(arr, rgb, tol=COLOR_TOL):
    r, g, b = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    return ((abs(r - rgb[0]) < tol) & (abs(g - rgb[1]) < tol) & (abs(b - rgb[2]) < tol))


def extract_png(pdf_path, page_num=4):
    """Pull the embedded Figure 2 PNG out of the PDF."""
    doc = fitz.open(pdf_path)
    page = doc[page_num]
    xref = page.get_images(full=True)[0][0]
    info = doc.extract_image(xref)
    return info["image"]


def measure_panel(arr, masks, centers, y_scan):
    """Mean segment height (in cases) for each bar in one panel."""
    blue_m, green_m, red_m = masks
    any_m = blue_m | green_m | red_m
    y0, y1 = y_scan
    out = []
    for cx in centers:
        lo, hi = int(round(cx - 15)), int(round(cx + 15))
        col_has = any_m[y0:y1, lo:hi].any(0)
        width = int(col_has.sum())
        if width == 0:
            out.append((0.0, 0.0, 0.0))
            continue

        def height(m):
            return m[y0:y1, lo:hi][:, col_has].sum() / width / PX_PER_UNIT

        out.append((height(blue_m), height(green_m), height(red_m)))
    return out


def main():
    root = Path(__file__).resolve().parent.parent
    png_bytes = extract_png(str(root / "papers" / "Andrejko 2024.pdf"))

    tmp_png = Path(__file__).resolve().parent / "_andrejko_fig2.png"
    tmp_png.write_bytes(png_bytes)
    arr = np.array(Image.open(tmp_png).convert("RGB")).astype(int)
    tmp_png.unlink()

    masks = (color_mask(arr, BLUE), color_mask(arr, GREEN), color_mask(arr, RED))

    n = len(SEROTYPES)
    spacing = (X_LAST - X_FIRST) / (n - 1)
    centers = [X_FIRST + spacing * k for k in range(n)]

    panels = [
        ("Early (May 2010-May 2014)", measure_panel(arr, masks, centers, PANEL_A["y_scan"])),
        ("Late (Jun 2014-Dec 2019)", measure_panel(arr, masks, centers, PANEL_B["y_scan"])),
    ]

    rows = []
    for period, data in panels:
        for serotype, (blue, green, red) in zip(SEROTYPES, data):
            d3, d12, d0 = round(blue), round(green), round(red)
            rows.append([period, serotype, "No PCV doses", d0, STUDY_TOTAL_NO_PCV])
            rows.append([period, serotype, "1-2 PCV13 doses", d12, STUDY_TOTAL_1_2_DOSES])
            rows.append([period, serotype, "3+ PCV13 doses", d3, STUDY_TOTAL_3PLUS_DOSES])

    out_path = root / "data" / "andrejko_figure_serotypes.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Study_Period", "Serotype", "Vaccine_Status", "Cases", "Total_Cases"])
        writer.writerows(rows)

    total = sum(r[3] for r in rows)
    s19a = sum(r[3] for r in rows if r[1] == "19A")
    print(f"Wrote {out_path} ({len(rows)} rows)")
    print(f"Total cases: {total} (Fig 1 reports 1161)")
    print(f"Serotype 19A total: {s19a} (paper: 96)")


if __name__ == "__main__":
    main()
