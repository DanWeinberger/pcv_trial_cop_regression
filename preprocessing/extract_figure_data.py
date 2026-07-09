"""
Extract data from the Whitney 2006 figure (Lancet 368:1497) using the pdf2plot
vector-geometry method (https://github.com/adamkucharski/pdf2plot).

Rather than rasterizing the page and doing pixel/image processing, this reads the
PDF's actual vector path coordinates: each bar is a true rectangle in the PDF
geometry. Bars are separated into two series by fill colour, the y-axis is
calibrated from the tick-label text positions, and pixel(point)->data conversion
is linear.

Chart structure: bars are OVERLAID (blue and green share each x-position), not
grouped. Blue = TOTAL cases; green (unvaccinated) is drawn in front. Therefore
"at least one dose" = blue - green.

Output is in long/tidy format: one row per (serotype, vaccine status), with
Total_Cases giving the study-wide denominator (cases + controls from Table 1)
for that vaccine status, for computing per-serotype rates.

Validation against the paper's text:
  - sum of totals ~= 776 (paper reports 782 enrolled cases)
  - sum of "at least one dose" over the 7 vaccine types = 65 (paper: exactly 65)
  - tallest bar is 19A ("the most common cause of disease")
  - most common vaccine types are 14 and 19F (matches text)

Requires: pymupdf   (pip install pymupdf)
"""

import csv
from pathlib import Path

import fitz  # PyMuPDF

# Fill colours of the two bar series (read directly from the vector geometry)
BLUE = (0.7405508756637573, 0.8695964217185974, 0.9264667630195618)   # "At least one dose"
GREEN = (0.0, 0.4420233368873596, 0.24319829046726227)                # "Unvaccinated"

# Y-axis calibration from tick-label text positions on the figure page
Y_BASELINE = 682.5   # page y of data value 0 (also the bar baseline)
Y_TOP = 550.0        # page y of data value 100
SCALE = 100.0 / (Y_BASELINE - Y_TOP)   # cases per PDF point

BAR_MAX_WIDTH = 6.0  # data bars are w=4.8; legend swatches are w=6.8 -> excluded

# Serotype labels with their x-axis centre positions, in plotting order
LABELS = [
    ("14", 292.8), ("19F", 299.5), ("6B", 308.1), ("18C", 314.9), ("4", 324.9),
    ("23F", 330.7), ("9V", 339.2),
    ("19A", 345.9), ("6A", 354.6), ("9N", 362.3), ("23B", 369.4), ("18B", 377.1),
    ("18F", 385.1), ("23A", 392.6),
    ("33F", 400.7), ("3", 410.5), ("22F", 416.2), ("7F", 425.1), ("15C", 431.7),
    ("12F", 439.6), ("38", 448.2), ("15B", 455.0), ("10A", 462.5), ("35B", 470.5),
    ("1", 480.6), ("11A", 486.0), ("35F", 494.0), ("7C", 502.7), ("8", 511.5),
    ("17F", 517.4), ("NT", 525.7), ("15A", 532.6), ("20", 541.5), ("33A", 548.1),
    ("34", 557.1),
]

VACCINE = {"14", "19F", "6B", "18C", "4", "23F", "9V"}

# Study-wide denominators from Table 1 (cases + controls), not derived from the
# figure: >=1 dose PCV = 393 cases + 1690 controls; unvaccinated = 389 + 822.
STUDY_TOTAL_VACCINATED = 2083
STUDY_TOTAL_UNVACCINATED = 1211


def color_close(a, b, tol=0.02):
    return a is not None and all(abs(x - y) < tol for x, y in zip(a, b))


def nearest_serotype(cx):
    return min(LABELS, key=lambda lbl: abs(lbl[1] - cx))[0]


def extract(pdf_path, page_num=2):
    doc = fitz.open(pdf_path)
    page = doc[page_num]

    blue_by_serotype = {s: 0.0 for s, _ in LABELS}
    green_by_serotype = {s: 0.0 for s, _ in LABELS}

    for drawing in page.get_drawings():
        fill = drawing.get("fill")
        for item in drawing["items"]:
            if item[0] != "re":
                continue
            rect = item[1]
            if rect.width > BAR_MAX_WIDTH:      # skip legend swatches
                continue
            height = Y_BASELINE - rect.y0        # bar height in PDF points
            if height <= 0.3:                    # skip zero-height noise
                continue
            cx = (rect.x0 + rect.x1) / 2
            serotype = nearest_serotype(cx)
            value = height * SCALE
            if color_close(fill, BLUE):
                blue_by_serotype[serotype] += value
            elif color_close(fill, GREEN):
                green_by_serotype[serotype] += value

    rows = []
    for serotype, _ in LABELS:
        dose = max(0, round(blue_by_serotype[serotype] - green_by_serotype[serotype]))
        unvax = round(green_by_serotype[serotype])
        rows.append((serotype, ">=1 dose", dose, STUDY_TOTAL_VACCINATED))
        rows.append((serotype, "Unvaccinated", unvax, STUDY_TOTAL_UNVACCINATED))
    return rows


def main():
    pdf_path = Path(__file__).resolve().parent.parent / "papers" / "Whitney 2006.pdf"
    rows = extract(str(pdf_path))

    out_path = Path(__file__).resolve().parent.parent / "data" / "whitney_figure_serotypes.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Serotype", "Vaccine_Status", "Cases", "Total_Cases"])
        writer.writerows(rows)

    total_cases = sum(r[2] for r in rows)
    vax_dose = sum(r[2] for r in rows
                   if r[1] == ">=1 dose" and r[0] in VACCINE)
    print(f"Wrote {out_path}")
    print(f"Sum of totals: {total_cases} (paper reports 782)")
    print(f">=1 dose among vaccine types: {vax_dose} (paper reports 65)")


if __name__ == "__main__":
    main()
