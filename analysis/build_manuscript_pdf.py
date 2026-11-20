#!/usr/bin/env python3
"""Render the editable manuscript Markdown source as a journal-style PDF."""

from __future__ import annotations

import re
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript" / "MMP1_collagen_OI_adaptive_unwinding_manuscript.md"
OUTPUT = ROOT / "manuscript" / "MMP1_collagen_OI_adaptive_unwinding_manuscript.pdf"
def register_fonts() -> None:
    candidates = [
        (
            Path("C:/Windows/Fonts/times.ttf"),
            Path("C:/Windows/Fonts/timesbd.ttf"),
            Path("C:/Windows/Fonts/timesi.ttf"),
            Path("C:/Windows/Fonts/timesbi.ttf"),
        ),
        (
            Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"),
            Path("/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"),
        ),
    ]
    font_paths = next((paths for paths in candidates if all(path.exists() for path in paths)), None)
    if font_paths is None:
        raise FileNotFoundError("Times New Roman font files were not found on this system")
    regular, bold, italic, bold_italic = font_paths
    pdfmetrics.registerFont(TTFont("TimesNR", regular))
    pdfmetrics.registerFont(TTFont("TimesNR-Bold", bold))
    pdfmetrics.registerFont(TTFont("TimesNR-Italic", italic))
    pdfmetrics.registerFont(TTFont("TimesNR-BoldItalic", bold_italic))
    pdfmetrics.registerFontFamily(
        "TimesNR",
        normal="TimesNR",
        bold="TimesNR-Bold",
        italic="TimesNR-Italic",
        boldItalic="TimesNR-BoldItalic",
    )


def inline_markup(text: str) -> str:
    text = escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`(.+?)`", r'<font name="Courier">\1</font>', text)
    text = re.sub(r"\*(.+?)\*", r"<i>\1</i>", text)
    return text


def styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "ManuscriptTitle",
            parent=base["Title"],
            fontName="TimesNR-Bold",
            fontSize=17,
            leading=20,
            alignment=TA_LEFT,
            spaceAfter=7 * mm,
            textColor=colors.black,
        ),
        "subtitle": ParagraphStyle(
            "ManuscriptSubtitle",
            parent=base["Normal"],
            fontName="TimesNR-Italic",
            fontSize=10.5,
            leading=13,
            alignment=TA_LEFT,
            spaceAfter=5 * mm,
        ),
        "h1": ParagraphStyle(
            "HeadingOne",
            parent=base["Heading1"],
            fontName="TimesNR-Bold",
            fontSize=13,
            leading=15,
            spaceBefore=5 * mm,
            spaceAfter=2.5 * mm,
            keepWithNext=True,
            textColor=colors.black,
        ),
        "h2": ParagraphStyle(
            "HeadingTwo",
            parent=base["Heading2"],
            fontName="TimesNR-Bold",
            fontSize=10.5,
            leading=13,
            spaceBefore=3.5 * mm,
            spaceAfter=1.8 * mm,
            keepWithNext=True,
            textColor=colors.black,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName="TimesNR",
            fontSize=9.5,
            leading=12.1,
            alignment=TA_JUSTIFY,
            firstLineIndent=5 * mm,
            spaceAfter=0.7 * mm,
            allowWidows=0,
            allowOrphans=0,
        ),
        "body_flush": ParagraphStyle(
            "BodyFlush",
            parent=base["BodyText"],
            fontName="TimesNR",
            fontSize=9.5,
            leading=12.1,
            alignment=TA_JUSTIFY,
            spaceAfter=2.4 * mm,
            allowWidows=0,
            allowOrphans=0,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=base["BodyText"],
            fontName="TimesNR",
            fontSize=8.3,
            leading=10.2,
            alignment=TA_LEFT,
            spaceBefore=1.3 * mm,
            spaceAfter=3 * mm,
        ),
        "table": ParagraphStyle(
            "TableText",
            parent=base["BodyText"],
            fontName="TimesNR",
            fontSize=7.8,
            leading=9.2,
            alignment=TA_CENTER,
        ),
        "table_header": ParagraphStyle(
            "TableHeader",
            parent=base["BodyText"],
            fontName="TimesNR-Bold",
            fontSize=7.8,
            leading=9.2,
            alignment=TA_CENTER,
        ),
    }


def figure_flowables(markdown_line: str, source_dir: Path) -> list:
    match = re.match(r"!\[[^\]]*\]\(([^)]+)\)", markdown_line)
    if not match:
        return []
    path = (source_dir / match.group(1)).resolve()
    image = Image(str(path))
    max_width = 170 * mm
    max_height = 198 * mm
    if path.name == "figure_structural_metrics_by_generation.png":
        max_width = 170 * mm
        max_height = 128 * mm
    elif path.name == "figure_generation_28_structures.png":
        max_width = 170 * mm
        max_height = 130 * mm
    elif path.name == "figure_residue_resolved_structural_changes.png":
        max_width = 170 * mm
        max_height = 130 * mm
    elif path.name == "figure_mmp1_conformational_reorganisation.png":
        max_width = 170 * mm
        max_height = 147 * mm
    scale = min(max_width / image.imageWidth, max_height / image.imageHeight)
    image.drawWidth = image.imageWidth * scale
    image.drawHeight = image.imageHeight * scale
    image.hAlign = "CENTER"
    space_before = 0 if path.name == "figure_mmp1_conformational_reorganisation.png" else 2 * mm
    return [Spacer(1, space_before), image]


def parse_table(lines: list[str], index: int, style_map: dict[str, ParagraphStyle]) -> tuple[Table, int]:
    rows: list[list[str]] = []
    while index < len(lines) and lines[index].strip().startswith("|"):
        cells = [cell.strip() for cell in lines[index].strip().strip("|").split("|")]
        if not all(re.fullmatch(r"[-: ]+", cell) for cell in cells):
            rows.append(cells)
        index += 1
    formatted = []
    for row_index, row in enumerate(rows):
        style = style_map["table_header"] if row_index == 0 else style_map["table"]
        formatted.append([Paragraph(inline_markup(cell), style) for cell in row])
    available = 170 * mm
    if len(formatted[0]) == 4 and "Catalytic-HPX" in rows[0][1]:
        widths = [available * value for value in [0.16, 0.32, 0.26, 0.26]]
    else:
        widths = [available * value for value in [0.18, 0.15, 0.15, 0.14, 0.18, 0.20]]
    table = Table(
        formatted,
        colWidths=widths[: len(formatted[0])],
        repeatRows=1,
        splitByRow=0,
        hAlign="CENTER",
    )
    table.setStyle(
        TableStyle(
            [
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
            ]
        )
    )
    return table, index


def footer(canvas, document) -> None:
    canvas.saveState()
    canvas.setFont("TimesNR", 8)
    canvas.setFillColor(colors.HexColor("#555555"))
    canvas.drawString(20 * mm, 12 * mm, "OI-associated collagen packing vulnerability")
    canvas.drawRightString(A4[0] - 20 * mm, 12 * mm, str(document.page))
    canvas.restoreState()


def build() -> None:
    register_fonts()
    style_map = styles()
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    story: list = []
    index = 0
    indent_body = False
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line.startswith("# "):
            story.append(Paragraph(inline_markup(line[2:]), style_map["title"]))
        elif line.startswith("## "):
            heading = line[3:]
            story.append(Paragraph(inline_markup(heading), style_map["h1"]))
            indent_body = heading not in {"Abstract", "References"}
        elif line.startswith("### "):
            story.append(Paragraph(inline_markup(line[4:]), style_map["h2"]))
        elif line.startswith("!["):
            items = figure_flowables(line, SOURCE.parent)
            caption_index = index + 1
            while caption_index < len(lines) and not lines[caption_index].strip():
                caption_index += 1
            if caption_index < len(lines):
                caption_line = lines[caption_index].strip()
                if caption_line.startswith("*") and caption_line.endswith("*") and caption_line[1:-1].startswith("Figure"):
                    items.append(Paragraph(inline_markup(caption_line[1:-1]), style_map["caption"]))
                    index = caption_index
            story.append(KeepTogether(items))
        elif line.startswith("|"):
            table, index = parse_table(lines, index, style_map)
            story.extend([table, Spacer(1, 1.5 * mm)])
            continue
        elif line.startswith("*") and line.endswith("*") and not line.startswith("**"):
            content = line[1:-1]
            if content.startswith("Figure") or content.startswith("Table"):
                story.append(Paragraph(inline_markup(content), style_map["caption"]))
            else:
                story.append(Paragraph(f"<i>{inline_markup(content)}</i>", style_map["subtitle"]))
        else:
            body_style = style_map["body"] if indent_body else style_map["body_flush"]
            story.append(Paragraph(inline_markup(line), body_style))
        index += 1

    document = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=A4,
        rightMargin=20 * mm,
        leftMargin=20 * mm,
        topMargin=18 * mm,
        bottomMargin=19 * mm,
        title="Osteogenesis imperfecta-associated glycine substitutions encode local triple-helix packing vulnerability relevant to MMP-1-mediated matrix remodelling",
        author="Anthony Nash; Chase Harms; Susanta K Sarkar",
        subject="Adaptive atomistic molecular dynamics manuscript",
    )
    document.build(story, onFirstPage=footer, onLaterPages=footer)
    print(OUTPUT)


if __name__ == "__main__":
    build()
