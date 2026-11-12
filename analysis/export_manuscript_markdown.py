#!/usr/bin/env python3
"""Export the verified manuscript DOCX into a readable Markdown source file."""

from __future__ import annotations

import shutil
from pathlib import Path

from docx import Document
from docx.table import Table
from docx.text.paragraph import Paragraph


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "output" / "docx" / "MMP1_collagen_OI_adaptive_unwinding_manuscript.docx"
OUT_DIR = ROOT / "output" / "manuscript"
OUT = OUT_DIR / "MMP1_collagen_OI_adaptive_unwinding_manuscript.md"
FIGURE_SOURCE = ROOT / "analysis_results" / "unwinding_report"
FIGURES = [
    "figure_adaptive_unwinding_progression.png",
    "figure_worker_variability_common_endpoint.png",
]


def iter_blocks(document: Document):
    """Yield paragraphs and tables in document order."""
    for child in document.element.body.iterchildren():
        if child.tag.endswith("}p"):
            yield Paragraph(child, document)
        elif child.tag.endswith("}tbl"):
            yield Table(child, document)


def escape_cell(text: str) -> str:
    return " ".join(text.split()).replace("|", "\\|")


def render_table(table: Table) -> list[str]:
    rows = [[escape_cell(cell.text) for cell in row.cells] for row in table.rows]
    if not rows:
        return []
    lines = ["| " + " | ".join(rows[0]) + " |"]
    lines.append("| " + " | ".join("---" for _ in rows[0]) + " |")
    lines.extend("| " + " | ".join(row) + " |" for row in rows[1:])
    return lines


def run_text(paragraph: Paragraph) -> str:
    parts: list[str] = []
    for run in paragraph.runs:
        text = run.text
        if not text:
            continue
        if run.bold:
            text = f"**{text}**"
        if run.italic:
            text = f"*{text}*"
        parts.append(text)
    return "".join(parts).strip()


def export() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    figures_out = OUT_DIR / "figures"
    figures_out.mkdir(parents=True, exist_ok=True)
    for name in FIGURES:
        shutil.copy2(FIGURE_SOURCE / name, figures_out / name)

    document = Document(DOCX)
    lines: list[str] = []
    figure_index = 0
    for block in iter_blocks(document):
        if isinstance(block, Table):
            lines.extend(render_table(block))
            lines.append("")
            continue

        paragraph = block
        has_drawing = bool(paragraph._p.xpath(".//w:drawing"))
        if has_drawing:
            if figure_index >= len(FIGURES):
                raise RuntimeError("More document figures were found than expected")
            name = FIGURES[figure_index]
            lines.extend([f"![Figure {figure_index + 1}](figures/{name})", ""])
            figure_index += 1
            continue

        text = run_text(paragraph)
        if not text:
            continue
        style = paragraph.style.name if paragraph.style else ""
        if style == "Title":
            lines.extend([f"# {text}", ""])
        elif style == "Subtitle":
            lines.extend([f"*{text}*", ""])
        elif style == "Heading 1":
            lines.extend([f"## {text}", ""])
        elif style == "Heading 2":
            lines.extend([f"### {text}", ""])
        elif style == "Caption":
            lines.extend([f"*{text}*", ""])
        elif text.startswith("**Keywords  **"):
            keywords = text.removeprefix("**Keywords  **").strip()
            lines.extend([f"**Keywords:** {keywords}", ""])
        else:
            lines.extend([text, ""])

    if figure_index != len(FIGURES):
        raise RuntimeError(f"Expected {len(FIGURES)} figures, exported {figure_index}")
    OUT.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    export()
