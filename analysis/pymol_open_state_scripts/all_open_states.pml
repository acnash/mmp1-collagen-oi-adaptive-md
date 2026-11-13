# Wild type and three generation-28 open-state mutants in one PyMOL session.
# Run from the repository root with:
# pymol analysis/pymol_open_state_scripts/all_open_states.pml
#
# Four scenes are created: Wild_type, G978S, G984C and G987R.

reinitialize

load analysis_results/generation_28_render_inputs/wild_type_generation_028_selected.pdb, wild_type
load analysis_results/generation_28_render_inputs/G978S_generation_028_selected.pdb, G978S
load analysis_results/generation_28_render_inputs/G984C_generation_028_selected.pdb, G984C
load analysis_results/generation_28_render_inputs/G987R_generation_028_selected.pdb, G987R

# Display one continuous copy of each complex; no periodic or crystallographic copies are generated.
remove solvent
hide everything, all

# MMP-1 and collagen cartoons.
show cartoon, chain A and resi 80-448
color gray80, chain A and resi 80-448
show cartoon, chain B+C+D
color purple, chain B+C+D

# Wild-type glycine sites.
select wild_type_sites, wild_type and chain B and resi 978+984+987
show sticks, wild_type_sites
set stick_radius, 0.50, wild_type_sites
color gray70, wild_type_sites and elem C
color blue, wild_type_sites and elem N
color red, wild_type_sites and elem O
color yellow, wild_type_sites and elem S
color white, wild_type_sites and elem H
label wild_type and chain B and resi 978 and name CA, "Gly978"
label wild_type and chain B and resi 984 and name CA, "Gly984"
label wild_type and chain B and resi 987 and name CA, "Gly987"

# G978S substitution.
select G978S_site, G978S and chain B and resi 978
show sticks, G978S_site
set stick_radius, 0.50, G978S_site
color gray70, G978S_site and elem C
color blue, G978S_site and elem N
color red, G978S_site and elem O
color yellow, G978S_site and elem S
color white, G978S_site and elem H
label G978S and chain B and resi 978 and name CA, "G978S"

# G984C substitution.
select G984C_site, G984C and chain B and resi 984
show sticks, G984C_site
set stick_radius, 0.50, G984C_site
color gray70, G984C_site and elem C
color blue, G984C_site and elem N
color red, G984C_site and elem O
color yellow, G984C_site and elem S
color white, G984C_site and elem H
label G984C and chain B and resi 984 and name CA, "G984C"

# G987R substitution.
select G987R_site, G987R and chain B and resi 987
show sticks, G987R_site
set stick_radius, 0.50, G987R_site
color gray70, G987R_site and elem C
color blue, G987R_site and elem N
color red, G987R_site and elem O
color yellow, G987R_site and elem S
color white, G987R_site and elem H
label G987R and chain B and resi 987 and name CA, "G987R"

# Shared rendering settings.
bg_color white
set orthoscopic, on
set ray_opaque_background, on
set ray_shadows, off
set antialias, 2
set specular, 0.18
set shininess, 24
set ambient, 0.35
set direct, 0.65
set cartoon_sampling, 10
set cartoon_loop_radius, 0.28
set cartoon_tube_radius, 0.38
set label_color, black
set label_outline_color, white
set label_size, 14
set label_font_id, 7

# Store a separately framed scene for each structure.
disable all
enable wild_type
orient wild_type and chain A+B+C+D
zoom wild_type and chain A+B+C+D, 4
scene Wild_type, store

disable all
enable G978S
orient G978S and chain A+B+C+D
zoom G978S and chain A+B+C+D, 4
scene G978S, store

disable all
enable G984C
orient G984C and chain A+B+C+D
zoom G984C and chain A+B+C+D, 4
scene G984C, store

disable all
enable G987R
orient G987R and chain A+B+C+D
zoom G987R and chain A+B+C+D, 4
scene G987R, store

scene Wild_type, recall
deselect

