# DeepPBS
## Geometry invariant Deep learning of Protein-DNA structures for Binding Specificity prediction

![alt text](https://github.com/timkartar/DeepPBS/blob/main/run/figs/Fig1_white.png?raw=true)


## Installation

### Install pythonic dependencies
Pythonic dependencies for DeepPBS are listed on deeppbs_linux.yml
We recommend installation via `conda` packagement tool.

If you do not have `conda` please conda installation instructions [Here](https://docs.anaconda.com/free/anaconda/install/index.html)

Once you have installed conda you can install all the python dependencies by running the following:

`conda env create -f deeppbs_linux.yml`

This will create a conda environment named `deeppbs` and install the dependencies there.
Note: The installation is tested on linux systems with cuda11.3 and cuda11.6, you may have to adjust Pyorch version number based on your system.
The project was developed on PyG2.0.1, although future versions of PyG are backwards compatible as of now, but we cannot guarantee stability on all versions.
For more information refer installation pages for [PyTorch](https://pytorch.org/get-started/locally/) and [PyG](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html)

### Install third party packages

The preprocessing scripts depend on 3DNA and Curves, we have provided the packages required in `dependencies/bin` and how to source them in `run/process/proc_source.sh`. 
However, please refer to `x3dna-v2.3-linux-64bit/x3dna-v2.3/license.txt` for fair usage of this version of 3DNA software.

### Install DeepPBS

```
git clone https://github.com/timkartar/DeepPBS
cd deepPBS
pip install .
```

## Quickstart

Pre-trained models are provided with the package.

Example pipeline for processing and predicting is as below:
```
cd deeppbs/run/process

// process and predict the structures present in `deeppbs/run/process/pdb` directory
// see process_config.json and pred_config.json for parameters and path details (note some
parameters are unused) and `proc_source.sh` for required environment setup
// outputs will be generated in `deeppbs/run/process/output` directory

./process_and_predict.sh

// run interpretation on an example pdb file present in `deeppbs/run/process/pdb` (provided it has 
// been processed to an npz file in `deeppbs/run/process/npz`) see `interpret_config.json` 
// output wil be generated in `deeppbs/run/plot_scripts/interpret_output/`

./vis_interpret.sh 5x6g
```

![output](https://github.com/timkartar/DeepPBS/assets/16060117/ff99b40a-432b-43ff-a2c5-b0bde06c2db5)
