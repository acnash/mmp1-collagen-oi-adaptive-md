#!/usr/bin/env python3
from __future__ import annotations

import math
import shutil
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/Users/anthony/Documents/MMP1 OI")
SRC = Path("/Users/anthony/Downloads")
MUT = ROOT / "generated_mutants"
OUT = MUT / "salted_150mM_NaCl"

AVOGADRO_PER_NM3_MOLAR = 0.602214076
TARGET_MOLAR = 0.150


@dataclass
class GroAtom:
    resnr: int
    resname: str
    atom: str
    nr: int
    rest: str
    xyz: tuple[float, float, float]
    vel: tuple[float, float, float] | None


def parse_gro(path: Path) -> tuple[str, list[GroAtom], str]:
    lines = path.read_text().splitlines()
    atoms: list[GroAtom] = []
    for line in lines[2:-1]:
        vals = [float(line[i:i + 8]) for i in range(20, len(line), 8) if line[i:i + 8].strip()]
        vel = tuple(vals[3:6]) if len(vals) >= 6 else None
        atoms.append(GroAtom(
            resnr=int(line[:5]),
            resname=line[5:10].strip(),
            atom=line[10:15].strip(),
            nr=int(line[15:20]),
            rest=line[20:],
            xyz=tuple(vals[:3]),
            vel=vel,
        ))
    return lines[0], atoms, lines[-1]


def format_atom(a: GroAtom, nr: int) -> str:
    vel = "" if a.vel is None else f"{a.vel[0]:8.4f}{a.vel[1]:8.4f}{a.vel[2]:8.4f}"
    return f"{a.resnr % 100000:5d}{a.resname:<5s}{a.atom:>5s}{nr % 100000:5d}{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}{vel}\n"


def water_residues(atoms: list[GroAtom]) -> list[tuple[int, list[int]]]:
    by_res: dict[int, list[int]] = {}
    for i, atom in enumerate(atoms):
        if atom.resname == "SOL":
            by_res.setdefault(atom.resnr, []).append(i)
    waters = []
    for resnr, idxs in sorted(by_res.items()):
        names = sorted(atoms[i].atom for i in idxs)
        if names == ["HW1", "HW2", "OW"]:
            waters.append((resnr, idxs))
    return waters


def add_salt_by_replacing_waters(gro_path: Path, out_path: Path, n_pairs: int) -> tuple[int, int, int]:
    title, atoms, box = parse_gro(gro_path)
    waters = water_residues(atoms)
    if len(waters) < 2 * n_pairs:
        raise RuntimeError(f"not enough waters in {gro_path}")

    chosen = waters[-2 * n_pairs:]
    na_waters = chosen[:n_pairs]
    cl_waters = chosen[n_pairs:]
    replace: dict[int, str] = {}
    for _, idxs in na_waters:
        for i in idxs:
            replace[i] = "NA"
    for _, idxs in cl_waters:
        for i in idxs:
            replace[i] = "CL"

    new_atoms: list[GroAtom] = []
    for i, atom in enumerate(atoms):
        ion = replace.get(i)
        if ion is None:
            new_atoms.append(atom)
            continue
        if atom.atom != "OW":
            continue
        new_atoms.append(GroAtom(atom.resnr, ion, ion, 0, "", atom.xyz, atom.vel))

    out_path.write_text(
        title + f" + {n_pairs} NaCl pairs (150 mM target)\n"
        + f"{len(new_atoms):5d}\n"
        + "".join(format_atom(atom, nr) for nr, atom in enumerate(new_atoms, 1))
        + box + "\n"
    )

    sol = sum(1 for atom in new_atoms if atom.resname == "SOL" and atom.atom == "OW")
    na = sum(1 for atom in new_atoms if atom.resname == "NA")
    cl = sum(1 for atom in new_atoms if atom.resname == "CL")
    return sol, na, cl


def update_top(source_top: Path, out_top: Path, sol: int, na: int, cl: int) -> None:
    lines = source_top.read_text().splitlines()
    out = []
    in_mol = False
    seen_na = False
    for line in lines:
        if line.strip() == "[ molecules ]":
            in_mol = True
            out.append(line)
            continue
        if in_mol:
            parts = line.split()
            if len(parts) >= 2 and not parts[0].startswith(";"):
                if parts[0] == "SOL":
                    out.append(f"SOL               {sol}")
                    continue
                if parts[0] == "NA":
                    out.append(f"NA                {na}")
                    seen_na = True
                    continue
                if parts[0] == "CL":
                    if not seen_na:
                        out.append(f"NA                {na}")
                        seen_na = True
                    out.append(f"CL                {cl}")
                    continue
        out.append(line)
    if in_mol and not seen_na:
        out.append(f"NA                {na}")
    out_top.write_text("\n".join(out) + "\n")


def box_pairs(gro_path: Path) -> tuple[float, int]:
    box = [float(x) for x in gro_path.read_text().splitlines()[-1].split()[:3]]
    vol = box[0] * box[1] * box[2]
    pairs = int(round(TARGET_MOLAR * vol * AVOGADRO_PER_NM3_MOLAR))
    return vol, pairs


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        shutil.copy2(src, dst / src.name)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    systems = [
        ("wild_type", SRC / "NPT_eq (1).gro", SRC / "system.top", SRC),
        ("collagen_G978S", MUT / "collagen_G978S" / "NPT_eq_collagen_G978S.gro", MUT / "collagen_G978S" / "system.top", MUT / "collagen_G978S"),
        ("collagen_G984C", MUT / "collagen_G984C" / "NPT_eq_collagen_G984C.gro", MUT / "collagen_G984C" / "system.top", MUT / "collagen_G984C"),
        ("collagen_G987R", MUT / "collagen_G987R" / "NPT_eq_collagen_G987R.gro", MUT / "collagen_G987R" / "system.top", MUT / "collagen_G987R"),
    ]
    for name, gro, top, include_dir in systems:
        vol, pairs = box_pairs(gro)
        sub = OUT / name
        sub.mkdir(exist_ok=True)
        out_gro = sub / f"NPT_eq_{name}_150mM_NaCl.gro"
        sol, na, cl = add_salt_by_replacing_waters(gro, out_gro, pairs)
        update_top(top, sub / "system.top", sol, na, cl)
        for include in ["forcefield.itp", "topol_MMP1_active.itp", "collagen.itp", "tip3p.itp", "ions.itp"]:
            copy_if_exists(include_dir / include, sub)
            copy_if_exists(SRC / include, sub)
        for mutant_itp in include_dir.glob("collagen_G*.itp"):
            copy_if_exists(mutant_itp, sub)
        (sub / "README.txt").write_text(
            f"{name}: target NaCl {TARGET_MOLAR:.3f} M, volume {vol:.6f} nm^3, "
            f"added {pairs} Na and {pairs} Cl by replacing {2 * pairs} waters. "
            f"Final SOL={sol}, NA={na}, CL={cl}.\n"
        )
        print(name, "pairs", pairs, "SOL", sol, "NA", na, "CL", cl, "atoms", out_gro)


if __name__ == "__main__":
    main()
