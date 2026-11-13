#!/usr/bin/env python3
"""Build one self-contained PyMOL Python script containing all four structures."""

from __future__ import annotations

import base64
import gzip
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
STRUCTURE_DIR = ROOT / "analysis_results" / "generation_28_render_inputs"
OUTPUT = HERE / "all_open_states_standalone.py"

INPUTS = {
    "wild_type": STRUCTURE_DIR / "wild_type_generation_028_selected.pdb",
    "G978S": STRUCTURE_DIR / "G978S_generation_028_selected.pdb",
    "G984C": STRUCTURE_DIR / "G984C_generation_028_selected.pdb",
    "G987R": STRUCTURE_DIR / "G987R_generation_028_selected.pdb",
}


def encoded_pdb(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("REMARK"):
        lines[0] = "REMARK Continuous selected open-state structure for PyMOL display"
    pdb_bytes = ("\n".join(lines) + "\n").encode("utf-8")
    return base64.b64encode(gzip.compress(pdb_bytes, compresslevel=9)).decode("ascii")


payload = {name: encoded_pdb(path) for name, path in INPUTS.items()}

source = '''#!/usr/bin/env python3
"""Self-contained PyMOL display of wild type and three open-state mutants.

In PyMOL, choose File > Run Script and select this file, or enter:
run all_open_states_standalone.py
"""

import base64
import gzip

from pymol import cmd


PDB_DATA = __PDB_DATA__

SITE_SELECTIONS = {
    "wild_type": "chain B and resi 978+984+987",
    "G978S": "chain B and resi 978",
    "G984C": "chain B and resi 984",
    "G987R": "chain B and resi 987",
}


def colour_sticks(selection):
    cmd.show("sticks", selection)
    cmd.set("stick_radius", 0.50, selection)
    cmd.color("gray70", f"{selection} and elem C")
    cmd.color("blue", f"{selection} and elem N")
    cmd.color("red", f"{selection} and elem O")
    cmd.color("yellow", f"{selection} and elem S")
    cmd.color("white", f"{selection} and elem H")


cmd.reinitialize()

for object_name, compressed_data in PDB_DATA.items():
    pdb_text = gzip.decompress(base64.b64decode(compressed_data)).decode("utf-8")
    cmd.read_pdbstr(pdb_text, object_name)

cmd.remove("solvent")
cmd.hide("everything", "all")

for object_name in PDB_DATA:
    # MMP-1 protein is residues 80-448. Residues 449-454 are metal ions;
    # four calcium atoms are named CA and must not be interpreted as C-alpha atoms.
    mmp1 = f"{object_name} and chain A and resi 80-448"
    collagen = f"{object_name} and chain B+C+D"
    cmd.show("cartoon", mmp1)
    cmd.color("gray80", mmp1)
    cmd.set("cartoon_color", "gray80", mmp1)
    cmd.show("cartoon", collagen)
    cmd.color("purple", collagen)
    cmd.set("cartoon_color", "purple", collagen)
    colour_sticks(f"{object_name} and ({SITE_SELECTIONS[object_name]})")

cmd.label("wild_type and chain B and resi 978 and name CA", '"Gly978"')
cmd.label("wild_type and chain B and resi 984 and name CA", '"Gly984"')
cmd.label("wild_type and chain B and resi 987 and name CA", '"Gly987"')
cmd.label("G978S and chain B and resi 978 and name CA", '"G978S"')
cmd.label("G984C and chain B and resi 984 and name CA", '"G984C"')
cmd.label("G987R and chain B and resi 987 and name CA", '"G987R"')

cmd.bg_color("white")
cmd.set("orthoscopic", 1)
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("specular", 0.18)
cmd.set("shininess", 24)
cmd.set("ambient", 0.35)
cmd.set("direct", 0.65)
cmd.set("cartoon_sampling", 10)
cmd.set("cartoon_loop_radius", 0.28)
cmd.set("cartoon_tube_radius", 0.38)
cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
cmd.set("label_size", 14)
cmd.set("label_font_id", 7)

scene_names = {
    "wild_type": "Wild_type",
    "G978S": "G978S",
    "G984C": "G984C",
    "G987R": "G987R",
}

for object_name, scene_name in scene_names.items():
    cmd.disable("all")
    cmd.enable(object_name)
    cmd.orient(f"{object_name} and chain A+B+C+D")
    cmd.zoom(f"{object_name} and chain A+B+C+D", buffer=4.0)
    cmd.scene(scene_name, "store")

cmd.scene("Wild_type", "recall")
cmd.deselect()
'''

source = source.replace("__PDB_DATA__", repr(payload))
OUTPUT.write_text(source, encoding="utf-8")
print(OUTPUT)
