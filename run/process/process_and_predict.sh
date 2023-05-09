#!/bin/bash
source proc_source.sh
python ../process_co_crystal.py input.txt process_config.json --no_pwm
rm *.pdb
rm *.par
rm *.pqr
rm *.r3d
rm *.dat
rm *.log
ls npz > predict_input.txt
python ../predict.py predict_input.txt ./output/ -c pred_config.json

