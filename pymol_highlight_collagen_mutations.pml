python
from pymol import cmd

models = [
    ("WT", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/wild_type/NPT_eq_wild_type_150mM_NaCl.gro"),
    ("G978S", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/collagen_G978S/NPT_eq_collagen_G978S_150mM_NaCl.gro"),
    ("G984C", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/collagen_G984C/NPT_eq_collagen_G984C_150mM_NaCl.gro"),
    ("G987R", "/Users/anthony/Documents/MMP1 OI/generated_mutants/salted_150mM_NaCl/collagen_G987R/NPT_eq_collagen_G987R_150mM_NaCl.gro"),
]

cmd.reinitialize()

for obj, filename in models:
    cmd.load(filename, obj)

cmd.hide("everything")
cmd.remove("resn SOL+NA+CL")

for obj, _ in models:
    cmd.show("cartoon", obj)
    cmd.color("gray70", obj)

cmd.select("WT_glycine_sites", "WT and resn GLY and resi 978+984+987")
cmd.select("G978S_mutation", "G978S and resn SER and resi 978")
cmd.select("G984C_mutation", "G984C and resn CYS and resi 984")
cmd.select("G987R_mutation", "G987R and resn ARG and resi 987")

cmd.show("sticks", "WT_glycine_sites or G978S_mutation or G984C_mutation or G987R_mutation")
cmd.show("spheres", "name CA and (WT_glycine_sites or G978S_mutation or G984C_mutation or G987R_mutation)")

cmd.color("yellow", "WT_glycine_sites")
cmd.color("cyan", "G978S_mutation")
cmd.color("orange", "G984C_mutation")
cmd.color("magenta", "G987R_mutation")

cmd.set("sphere_scale", 0.35, "WT_glycine_sites or G978S_mutation or G984C_mutation or G987R_mutation")
cmd.set("stick_radius", 0.16)
cmd.bg_color("white")

cmd.label("name CA and WT_glycine_sites", '"WT Gly" + resi')
cmd.label("name CA and G978S_mutation", '"G978S"')
cmd.label("name CA and G984C_mutation", '"G984C"')
cmd.label("name CA and G987R_mutation", '"G987R"')

cmd.orient("WT_glycine_sites or G978S_mutation or G984C_mutation or G987R_mutation")
cmd.zoom("WT_glycine_sites or G978S_mutation or G984C_mutation or G987R_mutation", 8)
python end

set label_size, 18
set ray_opaque_background, off
