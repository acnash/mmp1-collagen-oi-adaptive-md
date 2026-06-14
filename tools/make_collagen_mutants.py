#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import re
import shutil
from dataclasses import dataclass, replace
from pathlib import Path


ROOT = Path("/Users/anthony/Documents/MMP1 OI")
SRC = Path("/Users/anthony/Downloads")
OUT = ROOT / "generated_mutants"

COLLAGEN = SRC / "collagen.itp"
MMP = SRC / "topol_MMP1_active.itp"
GRO = SRC / "NPT_eq (1).gro"
TOP = SRC / "system.top"


@dataclass
class AtomRec:
    nr: int
    atype: str
    resnr: int
    resname: str
    atom: str
    cgnr: int
    charge: str
    mass: str
    comment: str = ""


@dataclass
class GroAtom:
    resnr: int
    resname: str
    atom: str
    nr: int
    xyz: tuple[float, float, float]
    vel: tuple[float, float, float] | None


def parse_atoms(lines: list[str]) -> tuple[list[str], list[AtomRec], list[str], dict[str, list[list[str]]]]:
    pre: list[str] = []
    atoms: list[AtomRec] = []
    sections: dict[str, list[list[str]]] = {}
    current = None
    in_atoms = False
    after_atoms = False
    for line in lines:
        m = re.match(r"\s*\[\s*([^\]]+)\s*\]", line)
        if m:
            name = m.group(1).strip()
            current = name
            if name == "atoms":
                in_atoms = True
                pre.append(line)
            else:
                if in_atoms:
                    in_atoms = False
                    after_atoms = True
                if after_atoms:
                    sections.setdefault(name, []).append([line])
                else:
                    pre.append(line)
            continue
        if in_atoms:
            if not line.strip() or line.lstrip().startswith(";"):
                pre.append(line)
                continue
            parts = line.split(";")
            fields = parts[0].split()
            if len(fields) >= 8 and fields[0].isdigit():
                atoms.append(AtomRec(
                    nr=int(fields[0]), atype=fields[1], resnr=int(fields[2]),
                    resname=fields[3], atom=fields[4], cgnr=int(fields[5]),
                    charge=fields[6], mass=fields[7],
                    comment=(";" + ";".join(parts[1:]).rstrip("\n")) if len(parts) > 1 else "",
                ))
            else:
                pre.append(line)
            continue
        if after_atoms and current:
            sections.setdefault(current, []).append([line])
        else:
            pre.append(line)
    return pre, atoms, [], sections


def flatten_section(section_lines: list[list[str]]) -> list[str]:
    return [x for sub in section_lines for x in sub]


def parse_numeric_section(section_lines: list[str]) -> tuple[list[str], list[list[str]]]:
    header: list[str] = []
    entries: list[list[str]] = []
    for line in section_lines:
        s = line.strip()
        if not s or s.startswith(";") or s.startswith("["):
            header.append(line)
            continue
        fields = s.split()
        if fields and fields[0].isdigit():
            entries.append(fields)
        else:
            header.append(line)
    return header, entries


def atom_line(a: AtomRec) -> str:
    return (
        f"{a.nr:6d} {a.atype:>11s} {a.resnr:7d} {a.resname:>6s} "
        f"{a.atom:>6s} {a.cgnr:6d} {a.charge:>10s} {a.mass:>10s}"
        + (f"     {a.comment}" if a.comment else "")
        + "\n"
    )


def section_entry(fields: list[str]) -> str:
    return " ".join(f"{x:>5s}" for x in fields) + "\n"


def residue_template(atoms: list[AtomRec], residue: str, resnr: int) -> list[AtomRec]:
    out = [a for a in atoms if a.resnr == resnr and a.resname == residue]
    if not out:
        raise RuntimeError(f"missing template {residue} {resnr}")
    return out


CYS_TEMPLATE = [
    ("N", "N", "-0.4157", "14.01"),
    ("H", "H", "0.2719", "1.008"),
    ("CT", "CA", "0.0213", "12.01"),
    ("H1", "HA", "0.1124", "1.008"),
    ("CT", "CB", "-0.1231", "12.01"),
    ("H1", "HB1", "0.1112", "1.008"),
    ("H1", "HB2", "0.1112", "1.008"),
    ("SH", "SG", "-0.3119", "32.06"),
    ("HS", "HG", "0.1933", "1.008"),
    ("C", "C", "0.5973", "12.01"),
    ("O", "O", "-0.5679", "16.0"),
]


def make_new_residue(template: list[AtomRec] | None, mut: str, resnr: int) -> list[AtomRec]:
    if mut == "CYS":
        return [AtomRec(0, t, resnr, "CYS", name, 0, chg, mass) for t, name, chg, mass in CYS_TEMPLATE]
    assert template is not None
    return [replace(a, nr=0, resnr=resnr, resname=mut, cgnr=0, comment="") for a in template]


def build_mapping(old_atoms: list[AtomRec], new_res: list[AtomRec], target_resnr: int) -> tuple[list[AtomRec], dict[int, int], dict[str, int]]:
    target_old = [a for a in old_atoms if a.resnr == target_resnr]
    old_by_name = {a.atom: a for a in target_old}
    new_by_name = {a.atom: a for a in new_res}
    old_to_new: dict[int, int] = {}
    new_atoms: list[AtomRec] = []
    next_nr = 1
    inserted = False
    for a in old_atoms:
        if a.resnr == target_resnr:
            if not inserted:
                for na in new_res:
                    na = replace(na, nr=next_nr, cgnr=next_nr)
                    new_atoms.append(na)
                    if na.atom in old_by_name:
                        old_to_new[old_by_name[na.atom].nr] = next_nr
                    elif na.atom == "HA" and "HA1" in old_by_name:
                        old_to_new[old_by_name["HA1"].nr] = next_nr
                    next_nr += 1
                inserted = True
            continue
        na = replace(a, nr=next_nr, cgnr=next_nr)
        new_atoms.append(na)
        old_to_new[a.nr] = next_nr
        next_nr += 1
    new_name_to_nr = {a.atom: a.nr for a in new_atoms if a.resnr == target_resnr}
    return new_atoms, old_to_new, new_name_to_nr


def atom_context(atoms: list[AtomRec], resnr: int) -> dict[tuple[str, str], int]:
    ctx: dict[tuple[str, str], int] = {}
    for a in atoms:
        role = "cur" if a.resnr == resnr else "prev" if a.resnr == resnr - 1 else "next" if a.resnr == resnr + 1 else ""
        if role:
            ctx[(role, a.atom)] = a.nr
    return ctx


def add_template_interactions(entries: list[list[str]], old_atoms: list[AtomRec], new_atoms: list[AtomRec],
                              donor_res: int, target_res: int, old_to_new: dict[int, int],
                              section_name: str) -> list[list[str]]:
    old_by_nr = {a.nr: a for a in old_atoms}
    donor_ctx = atom_context(old_atoms, donor_res)
    target_ctx = atom_context(new_atoms, target_res)
    donor_rev = {v: k for k, v in donor_ctx.items()}
    touched_roles = {"cur"}
    added: list[list[str]] = []
    for fields in entries:
        atom_count = 4 if section_name == "dihedrals" else 3 if section_name == "angles" else 2
        nums = [int(x) for x in fields[:atom_count]]
        if not any(n in donor_rev and donor_rev[n][0] in touched_roles for n in nums):
            continue
        mapped = []
        ok = True
        for n in nums:
            if n in donor_rev:
                key = donor_rev[n]
                if key in target_ctx:
                    mapped.append(target_ctx[key])
                else:
                    ok = False
                    break
            elif n in old_to_new:
                mapped.append(old_to_new[n])
            else:
                ok = False
                break
        if ok:
            added.append([str(x) for x in mapped] + fields[atom_count:])
    return added


def mutate_itp(mut_name: str, target_res: int, donor_res: int | None, donor_resname: str | None) -> tuple[str, int]:
    lines = COLLAGEN.read_text().splitlines(True)
    pre, atoms, _, sections_raw = parse_atoms(lines)
    mmp_pre, mmp_atoms, _, mmp_sections_raw = parse_atoms(MMP.read_text().splitlines(True))
    if donor_resname == "SER":
        tmpl = residue_template(mmp_atoms, "SER", 123)
    elif donor_resname == "ARG":
        tmpl = residue_template(atoms, "ARG", donor_res or 989)
    elif donor_resname == "CYS":
        tmpl = None
    else:
        raise RuntimeError(donor_resname)
    new_res = make_new_residue(tmpl, donor_resname, target_res)
    new_atoms, old_to_new, _ = build_mapping(atoms, new_res, target_res)
    old_target_nrs = {a.nr for a in atoms if a.resnr == target_res}
    old_name = "collagen"

    out = []
    for line in pre:
        out.append(line.replace(old_name, mut_name) if "[ moleculetype" not in line else line)
    out = [line.replace("collagen\t3", f"{mut_name}\t3") for line in out]
    out.extend(atom_line(a) for a in new_atoms)

    for section_name in ["bonds", "pairs", "angles", "dihedrals"]:
        raw = flatten_section(sections_raw.get(section_name, []))
        header, entries = parse_numeric_section(raw)
        out.extend(header)
        atom_count = 4 if section_name == "dihedrals" else 3 if section_name == "angles" else 2
        kept: list[list[str]] = []
        for fields in entries:
            nums = [int(x) for x in fields[:atom_count]]
            if any(n in old_target_nrs for n in nums):
                backbone_ok = True
                mapped = []
                for n in nums:
                    if n in old_target_nrs:
                        old_atom = next(a for a in atoms if a.nr == n)
                        if old_atom.atom == "HA2":
                            backbone_ok = False
                            break
                    if n in old_to_new:
                        mapped.append(old_to_new[n])
                    else:
                        backbone_ok = False
                        break
                if backbone_ok:
                    kept.append([str(x) for x in mapped] + fields[atom_count:])
                continue
            kept.append([str(old_to_new[int(x)]) for x in fields[:atom_count]] + fields[atom_count:])

        donor_entries = entries
        donor_atoms = atoms
        donor = donor_res
        if donor_resname == "SER":
            donor_raw = flatten_section(mmp_sections_raw.get(section_name, []))
            _, donor_entries = parse_numeric_section(donor_raw)
            donor_atoms = mmp_atoms
            donor = 123
        elif donor_resname == "CYS":
            donor_raw = flatten_section(mmp_sections_raw.get(section_name, []))
            _, donor_entries = parse_numeric_section(donor_raw)
            donor_atoms = mmp_atoms
            donor = 259
        added = add_template_interactions(donor_entries, donor_atoms, new_atoms, donor, target_res, old_to_new, section_name)
        if donor_resname == "CYS":
            names = {a.atom: a.nr for a in new_atoms if a.resnr == target_res}
            if section_name == "bonds":
                added.append([str(names["SG"]), str(names["HG"]), "1"])
            elif section_name == "pairs":
                for n in ["CA", "HB1", "HB2"]:
                    added.append([str(names[n]), str(names["HG"]), "1"])
            elif section_name == "angles":
                added.append([str(names["CB"]), str(names["SG"]), str(names["HG"]), "1"])
            elif section_name == "dihedrals":
                for a in ["CA", "HB1", "HB2"]:
                    added.append([str(names[a]), str(names["CB"]), str(names["SG"]), str(names["HG"]), "9"])
        seen = set()
        merged: list[list[str]] = []
        for e in kept + added:
            key = (section_name, tuple(e))
            if key not in seen:
                seen.add(key)
                merged.append(e)
        out.extend(section_entry(e) for e in merged)
    return "".join(out), len(new_atoms)


def parse_gro(path: Path) -> tuple[str, list[GroAtom], str]:
    lines = path.read_text().splitlines()
    title = lines[0]
    atoms: list[GroAtom] = []
    for line in lines[2:-1]:
        resnr = int(line[:5])
        resname = line[5:10].strip()
        atom = line[10:15].strip()
        nr = int(line[15:20])
        vals = [float(line[i:i+8]) for i in range(20, len(line), 8) if line[i:i+8].strip()]
        xyz = tuple(vals[:3])
        vel = tuple(vals[3:6]) if len(vals) >= 6 else None
        atoms.append(GroAtom(resnr, resname, atom, nr, xyz, vel))
    return title, atoms, lines[-1]


def format_gro_atom(a: GroAtom, nr: int) -> str:
    vel = "" if a.vel is None else f"{a.vel[0]:8.4f}{a.vel[1]:8.4f}{a.vel[2]:8.4f}"
    return f"{a.resnr % 100000:5d}{a.resname:<5s}{a.atom:>5s}{nr % 100000:5d}{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}{vel}\n"


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vadd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vmul(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n == 0:
        raise RuntimeError("zero-length vector")
    return vmul(a, 1.0 / n)


def frame(n_xyz, ca_xyz, c_xyz):
    e1 = unit(vsub(c_xyz, ca_xyz))
    nvec = vsub(n_xyz, ca_xyz)
    e3 = unit(cross(e1, nvec))
    e2 = cross(e3, e1)
    return ca_xyz, e1, e2, e3


def to_local(p, fr):
    origin, e1, e2, e3 = fr
    d = vsub(p, origin)
    return (dot(d, e1), dot(d, e2), dot(d, e3))


def from_local(local, fr):
    origin, e1, e2, e3 = fr
    return vadd(origin, vadd(vmul(e1, local[0]), vadd(vmul(e2, local[1]), vmul(e3, local[2]))))


def rotate_point(p, a, b, deg):
    theta = math.radians(deg)
    k = unit(vsub(b, a))
    v = vsub(p, a)
    term1 = vmul(v, math.cos(theta))
    term2 = vmul(cross(k, v), math.sin(theta))
    term3 = vmul(k, dot(k, v) * (1.0 - math.cos(theta)))
    return vadd(a, vadd(term1, vadd(term2, term3)))


def min_heavy_contact(res_atoms: list[GroAtom], env: list[GroAtom], names: set[str]) -> float:
    best = 999.0
    for a in res_atoms:
        if a.atom not in names or a.atom.startswith("H"):
            continue
        for b in env:
            if b.atom.startswith("H"):
                continue
            d = norm(vsub(a.xyz, b.xyz))
            if d < best:
                best = d
    return best


def optimize_sidechain(res_atoms: list[GroAtom], env: list[GroAtom], mut: str) -> list[GroAtom]:
    by = {a.atom: a for a in res_atoms}
    angles = [0, 60, 120, 180, 240, 300]
    if mut == "SER":
        stages = [(("CA", "CB"), {"OG", "HG"})]
        score_names = {"CB", "OG"}
    elif mut == "CYS":
        stages = [(("CA", "CB"), {"SG", "HG"})]
        score_names = {"CB", "SG"}
    elif mut == "ARG":
        stages = [
            (("CA", "CB"), {"CG", "HG1", "HG2", "CD", "HD1", "HD2", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"}),
            (("CB", "CG"), {"CD", "HD1", "HD2", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"}),
            (("CG", "CD"), {"NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"}),
        ]
        score_names = {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}
    else:
        return res_atoms

    best_atoms = res_atoms
    best_score = min_heavy_contact(res_atoms, env, score_names)

    def apply_stage(atoms_in, stage, angle):
        axis_names, moving = stage
        lookup = {a.atom: a for a in atoms_in}
        a = lookup[axis_names[0]].xyz
        b = lookup[axis_names[1]].xyz
        out = []
        for atom in atoms_in:
            if atom.atom in moving:
                out.append(replace(atom, xyz=rotate_point(atom.xyz, a, b, angle)))
            else:
                out.append(atom)
        return out

    def walk(stage_i, atoms_in):
        nonlocal best_atoms, best_score
        if stage_i == len(stages):
            score = min_heavy_contact(atoms_in, env, score_names)
            if score > best_score:
                best_score = score
                best_atoms = atoms_in
            return
        for angle in angles:
            walk(stage_i + 1, apply_stage(atoms_in, stages[stage_i], angle))

    walk(0, res_atoms)
    return best_atoms


def mutate_gro(target_res: int, mut: str, new_atom_names: list[str], arg_neutralize: bool) -> tuple[str, int, int]:
    title, atoms, box = parse_gro(GRO)
    mmp_n = 5825
    col_n = 443
    chain_start = mmp_n
    chain_end = mmp_n + col_n
    chain = atoms[chain_start:chain_end]
    donor_global_res = {"SER": 123, "ARG": 989, "CYS": 259}[mut]
    donor_pool = atoms if mut != "SER" and mut != "CYS" else atoms[:mmp_n]
    donor = [a for a in donor_pool if a.resnr == donor_global_res and a.resname == mut]
    if mut == "CYS" and not donor:
        donor = [a for a in atoms[:mmp_n] if a.resnr == 259 and a.resname == "CYS"]
    target_old = [a for a in chain if a.resnr == target_res]
    if not target_old:
        raise RuntimeError(f"target residue {target_res} not found")
    old_by = {a.atom: a for a in target_old}
    donor_by = {a.atom: a for a in donor}
    donor_frame = frame(donor_by["N"].xyz, donor_by["CA"].xyz, donor_by["C"].xyz)
    target_frame = frame(old_by["N"].xyz, old_by["CA"].xyz, old_by["C"].xyz)
    old_vel = old_by["CA"].vel
    new_res_atoms: list[GroAtom] = []
    for name in new_atom_names:
        if name in old_by and name != "HA":
            xyz = old_by[name].xyz
            vel = old_by[name].vel
        elif name == "HA" and "HA1" in old_by:
            xyz = old_by["HA1"].xyz
            vel = old_by["HA1"].vel
        elif name in donor_by:
            xyz = from_local(to_local(donor_by[name].xyz, donor_frame), target_frame)
            vel = old_vel
        elif mut == "CYS" and name == "HG":
            sg = next(x for x in new_res_atoms if x.atom == "SG").xyz
            cb = next(x for x in new_res_atoms if x.atom == "CB").xyz
            direction = unit(vsub(sg, cb))
            xyz = vadd(sg, vmul(direction, 0.134))
            vel = old_vel
        else:
            raise RuntimeError(f"no coordinate source for {mut} {name}")
        new_res_atoms.append(GroAtom(old_by["CA"].resnr, mut, name, 0, xyz, vel))

    env_atoms = [a for i, a in enumerate(atoms) if not (chain_start <= i < chain_end and a.resnr == target_res)]
    new_res_atoms = optimize_sidechain(new_res_atoms, env_atoms, mut)

    new_chain: list[GroAtom] = []
    inserted = False
    for a in chain:
        if a.resnr == target_res:
            if not inserted:
                new_chain.extend(new_res_atoms)
                inserted = True
            continue
        new_chain.append(a)

    new_atoms = atoms[:chain_start] + new_chain + atoms[chain_end:]
    sol_count = 18859
    cl_count = 10
    if arg_neutralize:
        # Replace the last water molecule before the chloride block with one chloride.
        sol_idxs = [i for i, a in enumerate(new_atoms) if a.resname == "SOL"]
        last_three = sol_idxs[-3:]
        ow = next(new_atoms[i] for i in last_three if new_atoms[i].atom == "OW")
        keep = [a for i, a in enumerate(new_atoms) if i not in set(last_three)]
        first_cl = next(i for i, a in enumerate(keep) if a.resname == "CL")
        added = GroAtom(ow.resnr, "CL", "CL", 0, ow.xyz, ow.vel)
        new_atoms = keep[:first_cl] + [added] + keep[first_cl:]
        sol_count -= 1
        cl_count += 1
    return write_gro_text(title + f" {target_res}{mut}", new_atoms, box), sol_count, cl_count


def write_gro_text(title: str, atoms: list[GroAtom], box: str) -> str:
    out = [title + "\n", f"{len(atoms):5d}\n"]
    for i, a in enumerate(atoms, 1):
        out.append(format_gro_atom(a, i))
    out.append(box + "\n")
    return "".join(out)


def write_top(mut_name: str, sol: int, cl: int) -> str:
    return f"""#include "forcefield.itp"

#include "topol_MMP1_active.itp"
;#include "posre_MMP1_active.itp"

#include "{mut_name}.itp"
#include "collagen.itp"
;#include "posre_collagen.itp"

#include "tip3p.itp"
#include "ions.itp"
[ system ]
Protein in water

[ molecules ]
MMP1_active       1
{mut_name:<16s} 1
collagen          2
SOL               {sol}
CL                {cl}
"""


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for dep in ["topol_MMP1_active.itp", "collagen.itp", "system.top"]:
        shutil.copy2(SRC / dep, OUT / f"original_{dep}")
    mutations = [
        ("collagen_G978S", 978, "SER", None, False),
        ("collagen_G984C", 984, "CYS", None, False),
        ("collagen_G987R", 987, "ARG", 989, True),
    ]
    for mut_name, res, residue, donor, neutralize in mutations:
        sub = OUT / mut_name
        sub.mkdir(exist_ok=True)
        shutil.copy2(MMP, sub / "topol_MMP1_active.itp")
        shutil.copy2(COLLAGEN, sub / "collagen.itp")
        itp, natoms = mutate_itp(mut_name, res, donor, residue)
        new_names = [a.atom for a in parse_atoms(itp.splitlines(True))[1] if a.resnr == res]
        gro, sol, cl = mutate_gro(res, residue, new_names, neutralize)
        (sub / f"{mut_name}.itp").write_text(itp)
        (sub / "system.top").write_text(write_top(mut_name, sol, cl))
        (sub / f"NPT_eq_{mut_name}.gro").write_text(gro)
        (sub / "README.txt").write_text(
            f"{mut_name}: first collagen chain only, residue {res} -> {residue}. "
            f"Mutant moleculetype atoms: {natoms}. SOL={sol}, CL={cl}.\n"
        )
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
