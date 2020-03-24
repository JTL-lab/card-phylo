# CARD-phylo

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.0.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/card-phylo.svg?branch=master)](https://travis-ci.org/snakemake-workflows/card-phylo)

This is a workflow to generate reference phylogenies for all AMR sequences in CARD canonical and prevalence.

Note the LICENSE file only refers to the workflow in this repository for licensing related to using CARD please visit the [website](https://card.mcmaster.ca/about).

Additionally if you make use of this please cite:

`Alcock et al. 2020. CARD 2020: Antibiotic Resistome Surveillance with the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 48, D517-D525.`

## Authors

* Finlay Maguire (@fmaguire)

## Usage

### Simple

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/card-phylo/releases).
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced). The latter way is recommended.

#### Step 2: Install and activate the conda environment

While you can set up individual envs for each rule this is a simple pipeline so we just have one big env.

`conda env create -f envs/card-phylo.yml`
`conda activate card-phylo`

#### Step 3: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

Main options are changing the version of CARD canonical/CARD prevalence and changing the cluster threshold proportion ID%.

#### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -cores $N -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

# Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/card-phylo) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.


## Testing

Tests cases are in the subfolder `.test`. They are automtically executed via continuous integration with Travis CI.
