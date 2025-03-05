[![DOI](https://img.shields.io/badge/DOI-NMETH-FFA500.svg)](https://doi.org/10.1038/s41592-024-02372-w)  

## Geometric deep learning of protein-DNA binding specificity

 <div>
  <a href="https://github.com/timkartar/DeepPBS">
         <img height="450" src="https://github.com/timkartar/DeepPBS/blob/main/run/figs/deeppbs.webp?raw=true" />
  </a>
 </div>

## Try it out on our web server (CPU only) 
 <div>
  <a href="https://rohslab.usc.edu/deeppbs/">
         <img height="50" src="https://static.vecteezy.com/system/resources/previews/021/351/649/original/web-server-icon-for-your-website-mobile-presentation-and-logo-design-free-vector.jpg" />
  </a>
 </div>
 
## Code ocean capsule
(same code structure as GitHub. If you wish to modify code/change input just copy the capsule into a new capsule that you own on Code Ocean!) 

[![DOI](https://8277274.fs1.hubspotusercontent-na1.net/hubfs/8277274/Code%20Ocean%20U4%20Theme%20Assets/code-ocean-footer-logo.svg)](https://doi.org/10.24433/CO.0545023.v2)

## Docker container
(requires a Linux machine with Docker installed. If you wish to use a GPU, please install the NVIDIA container toolkit as well)

run: `docker pull aricohen/deeppbs:latest`

wherever your terminal is, have a valid .cif or .pdb file available. Then, run the below command, replacing test.cif with your filename:
`docker run --gpus all -it -v $(pwd)/test.cif:/app/input/test.cif   -v $(pwd)/results:/output   deeppbs /app/input/test.cif`
(to only run prediction, and not generate heavy atom importance scores, add `-m` to the end of the above command)

This will create a folder named `results`, and place the important DeepPBS results inside. The `predict` folder contains results related to the prediction, including the position weight matrix. If you do not run -m, it will also generate an `interpretation` folder and put related results (Pymol session, residue wise scores) there.

## Installation
(should take 5-10 minutes with proper system setup)
### 1. Git clone the repository
```
git clone https://github.com/timkartar/DeepPBS
```
### 2. Install pythonic dependencies

We recommend installation via `conda` packagement tool.
If you do not have `conda` please refer conda installation instructions [Here](https://docs.anaconda.com/free/anaconda/install/index.html)

```
// gcc and cuda configs: gcc/12.3.0 cuda/12.2.1 (works with 12.2 and 12.1, just FYI)

conda create -n deeppbs_install python=3.10

conda init bash

conda activate deeppbs_install

// look here for other versions: https://pytorch.org/get-started/previous-versions/
conda install pytorch==2.3.0 torchvision==0.18.0 torchaudio==2.3.0 pytorch-cuda=12.1 -c pytorch -c nvidia

pip install torch_geometric

// look here for other versions: https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html
pip install torch_scatter torch_sparse torch_cluster -f https://data.pyg.org/whl/torch-2.3.0+cu121.html

pip install -U --no-cache-dir     biopython==1.83     logomaker     matplotlib==3.5.2     networkx     pandas==1.4.4     pdb2pqr     scipy==1.14.1     seaborn==0.13.2     freesasa==2.2.1 

```

### 3. Install DeepPBS

```
cd DeepPBS
pip install -e .
```
### 4. Third party packages

The preprocessing scripts depend on 3DNA and Curves, we have provided the packages required in `dependencies/bin` and how to source them in `run/process/proc_source.sh`. 
However, please refer to `x3dna-v2.3-linux-64bit/x3dna-v2.3/license.txt` for fair usage of this version of 3DNA software.

Note: The installation is tested on linux systems with cuda11.3 and cuda11.6, you may have to adjust Pyorch version number based on your system.
#### UPDATE (Feb 29, 2024): The latest version on github is tested on CUDA 12.2, PyTorch 2.3 and PyG 2.5. The `.yml` file has been updated accordingly.

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

-   `pymol` (opens pymol GUI)
-   `pip install matplotlib` (in the pymol GUI command prompt)
-    close the pymol GUI 
-   `pymol  ../plot_scripts/vis_interpret.py ../plot_scripts/ 5x6g.pdb` (run from terminal)

This will open a pymol session for the visualization (screenshot below) and save a .psw file in `run/plot_scripts/interpret_output`

![5x6g](https://github.com/timkartar/DeepPBS/blob/main/run/figs/5x6g.png?raw=true)

Simulation trajectories in PDB format snapshots can be processed in similar manner:

![output](https://github.com/timkartar/DeepPBS/blob/main/run/figs/output.gif?raw=true)

## Data availability

Figshare link: https://doi.org/10.6084/m9.figshare.25678053

## Run training

Download and place the data avilability number 2 somewhere on your system and configure the path in
`/run/config.json` (`"data_dir"`). Also configure the `"output_path"` as you wish.

run `./submit_cross.sh` . This will submit 5 cross-validation models to train simultaneaously.
Modify this script according to your need.
