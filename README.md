# This is the README.md for Daphnia_RestEggs_4lakes_smk_slurm

Snakemake workflow for the execution on the leo5 HPC cluster with the SLURM batch job sumission system. 


======================================================

## Information on the tutorial

This workflow is designated for the analysis of population structure, population summary statistics,  and gene flow, and the performance of ancestry paintings.

As input we use whole genome sequencing data in the form of realigned bam files. For more information on the pre-processing process see ...


======================================================


## Anaconda and Mamba


Mamba (https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++.

```
# Load Anaconda on cluster (here mach2):
module load Anaconda3/2023.10/miniconda-base-2023.10
 

# To use conda commands in your current shell session, first do:
eval "$(/$UIBK_CONDA_DIR/bin/conda shell.bash hook)"

# !!! NOT required for c770 group !!! 
## Create environment from yaml file (in envs/):
mamba env create -f envs/s21.yaml -p $SCRATCH/envs/daphnia
conda config --append envs_dirs $SCRATCH/envs

# Users of the c770 group can simply activate the environment on leo5 with:
conda activate daphnia

# Software can be installed/updated by modifying the envs/s21.yaml file.

# To use the modified conda environment, update with:
mamba env update --name daphnia --file envs/s21.yaml

```

## Get the snakemake pipeline for the slurm based cluster (here leo5)

```
mkdir <some_directory>
cd <some_directory>
rsync -avP /scratch/c7701178/snakemake/* ./

```


## Some information on the snakemake pipeline for the slurm based cluster (here leo5)

1. Testing and execution of snakemake in working directory (where your snakefile is located)

2. Slurm specific files and logs are in slurm/ directory
* cluster submission script for main job (slurm/clusterSnakemake.sh) **Do not forget to change the job-name and mail-user!**
* configuration profile for slurm (slurm/config.yaml)
* sbatch log files (output & error files) are sent to slurm/log/

3. The snakefile:
* contains the names of the output files (should correspond to the output in your rule)
* always includes the file rules/common.smk
* includes additional files with rules you want to execute (f.ex. rules/test.smk)

4. The rules/common.smk file:
* calls the input config file (config/config.yaml) and 
* contains python code for loading sample information (and helper functions)

5. The snakemake configuration file (config/config.yaml) contains the paths to adapters, Kraken2 database and references 

6. Sample information and metadata are in list/samples.tsv

7. Data was downloaded from https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz **NOT required for c770 group!**

## Local execution of snakemake (for testing only, else not recommended!)

```
# dry-run
snakemake -n

```

## Run snakemake via cluster execution on leo5

For more information on leo5 and slurm (https://www.uibk.ac.at/zid/systeme/hpc-systeme/common/tutorials/slurm-tutorial.html)

```
sbatch slurm/clusterSnakemake.sh

```

## Check jobs

```
# all jobs
sq
# your jobs only
squ
```

