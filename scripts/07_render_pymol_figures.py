#!/usr/bin/env python3
"""
Phase 5.2 - Render structural figures in PyMOL.

Run with: pymol -cq scripts/07_render_pymol_figures.py
(cq = no GUI, no splash; headless rendering)

Outputs (in figures/):
    pymol_full_length_ddg.png
    pymol_kinase_ddg.png
    pymol_kinase_conservation.png
    pymol_atp_allosteric.png
"""

from pymol import cmd
import os

FIG = "figures"
os.makedirs(FIG, exist_ok=True)


def reset_and_render(width=1400, height=1000):
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_transparency", 0.15)
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)


def render_full_length():
    cmd.reinitialize()
    reset_and_render()
    cmd.load("data/structures/alphafold_EGFR_with_ddg.pdb", "egfr")
    cmd.hide("everything")
    cmd.show("cartoon", "egfr")
    cmd.color("grey80", "egfr")
    # color by B-factor (ΔΔG + 3; range 0-15)
    cmd.spectrum("b", "blue_white_red", "egfr", minimum=0, maximum=15)
    # spheres at variant sites with significant ΔΔG (b > 3.5, i.e., ΔΔG > 0.5)
    cmd.show("spheres", "egfr and b > 6.0")
    cmd.set("sphere_scale", 0.4)
    cmd.orient()
    cmd.png(f"{FIG}/pymol_full_length_ddg.png",
             width=1400, height=1000, dpi=300, ray=1)
    print(f"wrote {FIG}/pymol_full_length_ddg.png")


def render_kinase_ddg():
    cmd.reinitialize()
    reset_and_render()
    cmd.load("data/structures/alphafold_EGFR_with_ddg.pdb", "egfr")
    cmd.hide("everything")
    cmd.show("cartoon", "egfr and resi 707-982")
    cmd.color("grey80", "egfr and resi 707-982")
    cmd.spectrum("b", "blue_white_red", "egfr and resi 707-982",
                 minimum=0, maximum=15)
    cmd.show("spheres", "egfr and resi 707-982 and b > 6.0")
    cmd.set("sphere_scale", 0.5)
    # highlight specific residues of interest
    cmd.show("sticks", "egfr and resi 790+858+719+861+768")
    cmd.color("black", "egfr and resi 790+858+719+861+768 and elem C")
    cmd.orient("egfr and resi 707-982")
    cmd.zoom("egfr and resi 707-982", buffer=5)
    cmd.png(f"{FIG}/pymol_kinase_ddg.png",
             width=1400, height=1000, dpi=300, ray=1)
    print(f"wrote {FIG}/pymol_kinase_ddg.png")

def render_kinase_conservation():
    """Discrete coloring: red = score 1, orange = score 2, green = score 3."""
    cmd.reinitialize()
    reset_and_render()
    cmd.load("data/structures/alphafold_EGFR_with_conservation.pdb", "egfr")
    cmd.hide("everything")
    cmd.show("cartoon", "egfr and resi 707-982")
    cmd.color("grey80", "egfr and resi 707-982")
    # discrete colors by conservation score, only on variant residues (b>0)
    cmd.color("red",    "(egfr) and (b > 0.0)")
    cmd.color("orange", "(egfr) and (b > 1.4)")
    cmd.color("green",  "(egfr) and (b > 2.4)")
    cmd.orient("egfr and resi 707-982")
    cmd.zoom("egfr and resi 707-982", buffer=5)
    cmd.png(f"{FIG}/pymol_kinase_conservation.png",
             width=1400, height=1000, dpi=300, ray=1)
    print(f"wrote {FIG}/pymol_kinase_conservation.png")


def render_atp_allosteric():
    cmd.reinitialize()
    reset_and_render()
    cmd.load("data/structures/alphafold_EGFR_with_ddg.pdb", "egfr")
    cmd.hide("everything")
    cmd.show("cartoon", "egfr and resi 707-982")
    cmd.color("grey70", "egfr and resi 707-982")
    # ATP site residues: K745, T790, M793, C797, D855, F856
    cmd.show("sticks", "egfr and resi 745+790+793+797+855+856")
    cmd.color("orange", "egfr and resi 745+790+793+797+855+856 and elem C")
    # allosteric pocket: L858
    cmd.show("sticks", "egfr and resi 858")
    cmd.color("purple", "egfr and resi 858 and elem C")
    cmd.orient("egfr and resi 745+790+793+797+855+856+858")
    cmd.zoom("egfr and resi 745+790+793+797+855+856+858", buffer=10)
    cmd.png(f"{FIG}/pymol_atp_allosteric.png",
             width=1400, height=1000, dpi=300, ray=1)
    print(f"wrote {FIG}/pymol_atp_allosteric.png")


def main():
    print("[Phase 5.2] Rendering PyMOL figures\n")
    render_full_length()
    render_kinase_ddg()
    render_kinase_conservation()
    render_atp_allosteric()
    print("\n[Phase 5.2] Done.")


main()
