python
from pymol import cmd

models = [
    ("WT", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/wild_type/NPT_eq_wild_type_150mM_NaCl.gro"),
    ("G978S", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/collagen_G978S/NPT_eq_collagen_G978S_150mM_NaCl.gro"),
    ("G984C", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/collagen_G984C/NPT_eq_collagen_G984C_150mM_NaCl.gro"),
    ("G987R", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/collagen_G987R/NPT_eq_collagen_G987R_150mM_NaCl.gro"),
]

cmd.reinitialize()

for obj, path in models:
    cmd.load(path, obj)

cmd.hide("everything")
cmd.remove("resn SOL")

for obj, _ in models:
    cmd.show("cartoon", obj)
    cmd.show("sticks", f"{obj} and not resn NA+CL+ZN+CA+M1+M2+M3+M4+M5+M6+SOL")
    cmd.color("gray75", obj)

cmd.select("all_ions", "resn NA+CL+ZN+CA+M1+M2+M3+M4+M5+M6")
cmd.show("spheres", "all_ions")
cmd.color("blue", "resn NA")
cmd.color("green", "resn CL")
cmd.color("slate", "resn CA+M3+M4+M5+M6")
cmd.color("orange", "resn ZN+M1+M2")

cmd.set("sphere_scale", 0.35, "resn NA+CL")
cmd.set("sphere_scale", 0.55, "resn ZN+CA+M1+M2+M3+M4+M5+M6")
cmd.set("cartoon_transparency", 0.15)
cmd.set("stick_radius", 0.08)
cmd.bg_color("white")

cmd.orient("not resn SOL")
cmd.zoom("not resn SOL", 4)
python end
