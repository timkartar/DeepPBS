#!/bin/bash
source proc_source.sh
python ../process_co_crystal.py input.txt process_config.json --no_pwm
rm *.pdb
rm *.par
rm *.pqr
rm *.r3d
rm *.dat
rm *.log
cat input.txt | sed 's/pdb/npz/g' | sed 's/cif/npz/g' > predict_input.txt
python ../predict.py predict_input.txt ./output/ -c ./pred_configs/pred_config_deeppbs.json

