# DeepPBS
## Geometric Deep learning on Protein-DNA structures for interpretable prediction of binding specificity across protein families

![alt text](https://github.com/timkartar/DeepPBS/blob/main/run/figs/Fig1_square_white.png?raw=true)


## Installation
(should take 5-10 minutes with proper system setup)
### 1. Git clone the repository
```
git clone https://github.com/timkartar/DeepPBS
```
### 2. Install pythonic dependencies

Pythonic dependencies for DeepPBS are listed on `deeppbs_linux.yml`. We recommend installation via `conda` packagement tool.
If you do not have `conda` please refer conda installation instructions [Here](https://docs.anaconda.com/free/anaconda/install/index.html)
```
cd DeepPBS
conda env create -f deeppbs_linux.yml
conda activate deeppbs
```
### 3. Install DeepPBS

```
pip install .
```
### 4. Third party packages

The preprocessing scripts depend on 3DNA and Curves, we have provided the packages required in `dependencies/bin` and how to source them in `run/process/proc_source.sh`. 
However, please refer to `x3dna-v2.3-linux-64bit/x3dna-v2.3/license.txt` for fair usage of this version of 3DNA software.

Note: The installation is tested on linux systems with cuda11.3 and cuda11.6, you may have to adjust Pyorch version number based on your system.
The project was developed on PyG2.0.1, although future versions of PyG are backwards compatible as of now, but we cannot guarantee stability on all versions.
For more information refer installation pages for [PyTorch](https://pytorch.org/get-started/locally/) and [PyG](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html)

## Usage pipeline for pre-trained DeepPBS

Example pipeline for processing and predicting is as below:

1. `cd run/process/`
2. Put your PDB files containing biological aseemblies of interest into `pdb` directory
3. run `ls pdb > input.txt`
4. `./process_and_predict.sh` (you can parallelize the steps in this script through multiple job submissions)

This will process the list of pdbs and put the processed npz files into `npz` directory.

Note: As evident, you can parallelize this script, but in that case make sure you create a separate working directory for each job. Otherwise temporary files generated during processing may conflict.

Then it will make predictions using the DeepPBS ensemble and put the predictions in `output` directory (in `run/process`)
Combined pre-processing and inference time for one biological assembly is in the order of seconds (e.g., for PDB ID 5x6g, about 15-20 seconds)

## Compute and Visualize perturbation based heavy atom interpretability
1. `cd run/process`
2. `./vis_interpret.sh <pdb_name_without .pdb>`, for example `./vis_interpret.sh 5x6g` 

This will compute and store the perturbation outcomes and other required information in `run/plot_scripts/interpret_output`

3.  You need a [PyMol](https://pymol.org/2/) executable for this step! Once installed, you can run the following

4.  `pymol  ../plot_scripts/vis_interpret.py ../plot_scripts/ 5x6g.pdb` 

This will open a pymol session for the visualization (screenshot below) and save a .psw file in `run/plot_scripts/interpret_output`

![5x6g](https://github.com/timkartar/DeepPBS/blob/main/run/figs/5x6g.png?raw=true)

Simulation trajectories in PDB format snapshots can be processed in similar manner:

![output](https://github.com/timkartar/DeepPBS/blob/main/run/figs/output.gif?raw=true)

## Data availability

1. The full set of PDB chain ID (which have DNA in the biological assembly) and corresponding PWM ids are available in a clusterwise manner in `data/jaspar_h11mo_cluster_wise_dna_containing_dataset.npy` (note: all of these do not pass processing criteria)
2. Processed npz files which may go as input to DeepPBS: [Here](https://drive.google.com/file/d/17kb7RgDFoD3nVBgjd1LELHCXTWX1bHnN/view?usp=sharing)
3. All analyzed Exd-Scr simulations frames in PDB format: [Here](https://drive.google.com/file/d/1-gcE0ykb-Bg_birUSPLY1vYOi138cgit/view?usp=sharing)
4. Cross-validation set is listed as `run/fold/valid*.txt` and benchmark set is listed as `run/folds/id.txt`
5. Simulation parameters: [Here](https://doi.org/10.6084/m9.figshare.22695778.v3)
6. List of protein sequence IDs used for application on predicted structures: [Here](https://github.com/timkartar/DeepPBS/tree/main/data/ids_for_predicted_structures)
7. Custom scripts and intermediate files used for data gathering from JASPAR, HOCOMOCO and PDB: [Here](https://drive.google.com/file/d/17kb7RgDFoD3nVBgjd1LELHCXTWX1bHnN/view?usp=drive_link)

## Acknowledgement
Parts of pre-processing code was contributed by Jared Sagendorf from previous projects in the Rohs Lab: DNAProDB (https://github.com/jaredsagendorf/dnaprodb-back) and GEMENAI (https://github.com/jaredsagendorf/gemenai)
