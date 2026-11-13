# G987R MMP-1/model-collagen complex: generation-28 selected open structure.
# Run from the repository root with: pymol analysis/pymol_open_state_scripts/G987R_open_state.pml

reinitialize
load analysis_results/generation_28_render_inputs/G987R_generation_028_selected.pdb, G987R

# Display one continuous complex only; no crystallographic or periodic copies are generated.
remove solvent
hide everything, all
select mmp1, G987R and chain A and resi 80-448
select collagen, G987R and chain B+C+D
select highlighted_residue, G987R and chain B and resi 987

show cartoon, mmp1 or collagen
color gray80, mmp1
set cartoon_color, gray80, mmp1
color purple, collagen
set cartoon_color, purple, collagen

show sticks, highlighted_residue
set stick_radius, 0.50, highlighted_residue
color gray70, highlighted_residue and elem C
color blue, highlighted_residue and elem N
color red, highlighted_residue and elem O
color yellow, highlighted_residue and elem S
color white, highlighted_residue and elem H
label G987R and chain B and resi 987 and name CA, "G987R"

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

orient G987R and chain A+B+C+D
zoom G987R and chain A+B+C+D, 4
deselect

