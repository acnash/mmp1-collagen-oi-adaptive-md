#!/usr/bin/env python3
"""Render aligned generation-28 MMP-1/collagen structures with PyMOL."""

import os
from pathlib import Path

from pymol import cmd


ROOT = Path.cwd().resolve()
STRUCTURE_DIR = ROOT / "analysis_results" / "generation_28_render_inputs"
OUTPUT_DIR = ROOT / "analysis_results" / "generation_28_pymol"
SYSTEMS = ["wild_type", "G978S", "G984C", "G987R"]
MUTATION_POSITION = {"G978S": 978, "G984C": 984, "G987R": 987}
WILD_TYPE_GLYCINES = [978, 984, 987]


def prepare_object(name: str) -> None:
    cmd.hide("everything", name)
    cmd.remove(f"{name} and hydro")

    enzyme = f"{name} and chain A"
    collagen = f"{name} and chain B+C+D"
    cmd.show("cartoon", enzyme)
    cmd.color("gray40", enzyme)
    cmd.set("cartoon_color", "gray40", enzyme)

    cmd.show("cartoon", collagen)
    for chain, colour in (("B", "teal"), ("C", "marine"), ("D", "yelloworange")):
        selection = f"{name} and chain {chain}"
        cmd.color(colour, selection)
        cmd.set("cartoon_color", colour, selection)

    if name == "wild_type":
        glycines = f"{name} and chain B and resi {'+'.join(map(str, WILD_TYPE_GLYCINES))}"
        cmd.show("sticks", glycines)
        cmd.color("orange", glycines)
        cmd.set("stick_color", "orange", glycines)
        cmd.set("stick_radius", 0.38, glycines)
        for position in WILD_TYPE_GLYCINES:
            cmd.label(f"{name} and chain B and resi {position} and name CA", f'"Gly{position}"')

    if name in MUTATION_POSITION:
        mutation = f"{name} and chain B and resi {MUTATION_POSITION[name]}"
        cmd.show("sticks", mutation)
        cmd.color("magenta", mutation)
        cmd.set("stick_color", "magenta", mutation)
        cmd.set("stick_radius", 0.42, mutation)
        cmd.label(f"{mutation} and name CA", f'"{name}"')


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cmd.reinitialize()
    for system in SYSTEMS:
        path = STRUCTURE_DIR / f"{system}_generation_028_selected.pdb"
        if not path.exists():
            raise FileNotFoundError(path)
        cmd.load(str(path), system)

    for system in SYSTEMS[1:]:
        cmd.super(f"{system} and chain A and name CA", "wild_type and chain A and name CA", cycles=5)

    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.18)
    cmd.set("shininess", 24)
    cmd.set("ambient", 0.34)
    cmd.set("direct", 0.66)
    cmd.set("cartoon_sampling", 10)
    cmd.set("cartoon_loop_radius", 0.28)
    cmd.set("cartoon_tube_radius", 0.38)
    cmd.set("label_color", "black")
    cmd.set("label_size", 22)
    cmd.set("label_font_id", 7)
    cmd.set("label_outline_color", "white")

    for system in SYSTEMS:
        prepare_object(system)
    requested = os.environ.get("PYMOL_SYSTEM")
    render_systems = [requested] if requested else SYSTEMS
    for system in render_systems:
        cmd.disable("all")
        cmd.enable(system)
        if system == "wild_type":
            focus = f"{system} and chain B+C+D and resi 976-989"
        else:
            position = MUTATION_POSITION[system]
            focus = f"{system} and chain B+C+D and resi {position - 2}-{position + 2}"
        cmd.orient(focus)
        cmd.zoom(focus, buffer=3.0)
        cmd.turn("y", -10)
        cmd.turn("x", 4)
        cmd.png(str(OUTPUT_DIR / f"{system}.png"), width=1500, height=1200, dpi=300, ray=0)

    if not requested:
        cmd.enable("all")
        cmd.save(str(OUTPUT_DIR / "generation_28_mmp1_collagen.pse"))
    cmd.quit()


main()
