# DeepPBS


Developement Environment:

* Biopython 1.79
* datetime 
* matplotlib 3.5.2
* scipy 1.7.3
* scikit-learn 1.1.1
* pytorch 1.12.1+cu116
* pytorch_geometric 2.1.0
* tqdm 4.64.0
* glob2 0.7
* igl 2.2.1
* json
* logomaker 0.8
* networkx 2.8.4
* pandas 1.4.4
* seaborn 0.11.2
* trimesh 3.10.0
* numpy 1.21.5

# Welcome to Deep DNAshape

The package includes an executable `deepDNAshape` to predict DNA shape features for any sequences.

It also includes all the components of Deep DNAshape. You may incoorporate Deep DNAshape into your pipeline or modify it to fit your needs.

## Installation

Prerequsite: `tensorflow >= 2.0` `numpy`
### Download and install through pip
```
git clone https://github.com/timkartar/DeepPBS
cd deepPBS
pip install .
```

## Quickstart

Pre-trained models are provided with the package.

Example pipeline for processing and predicting is as below:
```
cd deeppbs/bin/process

// process and predict the structures present in `deeppbs/bin/process/pdb` directory
// see process_config.json and pred_config.json for parameters and path details (note some
parameters are unused)
./process_and_predict.sh
// outputs will be generated in `deeppbs/bin/process/output` directory

// run interpretation on an example pdb file present in `deeppbs/bin/process/pdb` (provided it has 
// been processed to an npz file in `deeppbs/bin/process/npz`) see `interpret_config.json` 
./vis_interpret.sh 5x6g
// output wil be generated in `deeppbs/bin/plot_scripts/interpret_output/`

```
