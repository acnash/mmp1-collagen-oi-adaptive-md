# Wild-type MMP-1/model-collagen complex: generation-28 selected structure.
# Run from the repository root with: pymol analysis/pymol_open_state_scripts/wild_type_open_state.pml

reinitialize
load analysis_results/generation_28_render_inputs/wild_type_generation_028_selected.pdb, wild_type

# Display one continuous complex only; no crystallographic or periodic copies are generated.
remove solvent
hide everything, all
select mmp1, wild_type and chain A and resi 80-448
select collagen, wild_type and chain B+C+D
select highlighted_residues, wild_type and chain B and resi 978+984+987

show cartoon, mmp1 or collagen
color gray80, mmp1
set cartoon_color, gray80, mmp1
color purple, collagen
set cartoon_color, purple, collagen

show sticks, highlighted_residues
set stick_radius, 0.50, highlighted_residues
color gray70, highlighted_residues and elem C
color blue, highlighted_residues and elem N
color red, highlighted_residues and elem O
color yellow, highlighted_residues and elem S
color white, highlighted_residues and elem H

label wild_type and chain B and resi 978 and name CA, "Gly978"
label wild_type and chain B and resi 984 and name CA, "Gly984"
label wild_type and chain B and resi 987 and name CA, "Gly987"

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

orient wild_type and chain A+B+C+D
zoom wild_type and chain A+B+C+D, 4
deselect

