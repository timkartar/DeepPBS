#!/bin/bash
python ../interpret.py $1.npz -c interpret_config.json
/home/raktim/Downloads/pymol/bin/pymol  ../plot_scripts/vis_interpret.py ../plot_scripts/ $1.pdb
