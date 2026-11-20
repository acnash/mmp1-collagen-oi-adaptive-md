#!/usr/bin/env python3
"""Render the system-preparation Supplementary Information as a PDF."""

from __future__ import annotations

import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.platypus import Image, Paragraph, Spacer, Table, TableStyle

import build_manuscript_pdf as base


ROOT = Path(__file__).resolve().parents[1]
base.SOURCE = ROOT / "manuscript" / "MMP1_collagen_OI_adaptive_unwinding_supplementary_information.md"
base.OUTPUT = ROOT / "manuscript" / "MMP1_collagen_OI_adaptive_unwinding_supplementary_information.pdf"


def figure_flowables(markdown_line: str, source_dir: Path) -> list:
    match = re.match(r"!\[[^\]]*\]\(([^)]+)\)", markdown_line)
    if not match:
        return []
    path = (source_dir / match.group(1)).resolve()
    image = Image(str(path))
    height_by_name = {
        "figure_S1_preparation_workflow.png": 92 * mm,
        "figure_S2_local_cluster_minimisation.png": 105 * mm,
        "figure_S3_G978S_validation_thermodynamics.png": 137 * mm,
        "figure_S4_preparation_structural_change.png": 104 * mm,
        "figure_S5_final_state_quantities.png": 88 * mm,
    }
    max_width = 170 * mm
    max_height = height_by_name.get(path.name, 150 * mm)
    scale = min(max_width / image.imageWidth, max_height / image.imageHeight)
    image.drawWidth = image.imageWidth * scale
    image.drawHeight = image.imageHeight * scale
    image.hAlign = "CENTER"
    return [Spacer(1, 2 * mm), image]


def parse_table(lines: list[str], index: int, style_map) -> tuple[Table, int]:
    rows = []
    while index < len(lines) and lines[index].strip().startswith("|"):
        cells = [cell.strip() for cell in lines[index].strip().strip("|").split("|")]
        if not all(re.fullmatch(r"[-: ]+", cell) for cell in cells):
            rows.append(cells)
        index += 1
    formatted = []
    for row_index, row in enumerate(rows):
        style = style_map["table_header"] if row_index == 0 else style_map["table"]
        formatted.append([Paragraph(base.inline_markup(cell), style) for cell in row])
    ncols = len(formatted[0])
    available = 170 * mm
    first = 0.21 if ncols >= 4 else 0.28
    widths = [available * first] + [available * (1 - first) / (ncols - 1)] * (ncols - 1)
    table = Table(formatted, colWidths=widths, repeatRows=1, splitByRow=0, hAlign="CENTER")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E6E6E6")),
        ("TEXTCOLOR", (0, 0), (-1, -1), colors.black),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("ALIGN", (0, 1), (0, -1), "LEFT"),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 0.6, colors.black),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F6F6F6")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return table, index


def footer(canvas, document) -> None:
    canvas.saveState()
    canvas.setFont("TimesNR", 8)
    canvas.setFillColor(colors.HexColor("#555555"))
    canvas.drawString(20 * mm, 12 * mm, "Supplementary Information | OI-associated collagen packing vulnerability")
    canvas.drawRightString(A4[0] - 20 * mm, 12 * mm, f"S{document.page}")
    canvas.restoreState()


base.figure_flowables = figure_flowables
base.parse_table = parse_table
base.footer = footer


if __name__ == "__main__":
    base.build()
