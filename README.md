# CARD-phylo

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.11.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/fmaguire/card-phylo.svg?branch=master)](https://travis-ci.org/fmaguire/card-phylo)

This is a simplified version of the workflow to generate reference phylogenies for all AMR sequences in CARD canonical 

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

#### Step 2: Install dependencies

To run this you must have a newish [snakemake](https://snakemake.readthedocs.io/en/stable/) (>5.11) and [conda](https://docs.conda.io/en/latest/miniconda.html) available in your path.

This workflow has only been tested on linux systems, it should work on OSX and hopefully windows but no guarantees.

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config/config.yaml`.

Main options are changing the version of CARD canonical/CARD prevalence and changing the cluster threshold proportion ID%.

#### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -cores $N -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

#### Output 

    card-phylo_canonical_3.0.9/
    ├── card
    │   └──── canonical
    │       └── 3.0.9
    │           ├── card.json
    │           ├── protein_fasta_protein_homolog_model.fasta
    │           ├── protein_fasta_protein_knockout_model.fasta
    │           ├── protein_fasta_protein_overexpression_model.fasta
    │           └── protein_fasta_protein_variant_model.fasta
    ├── card_protein.fasta
    ├── mmseqs
    │   ├── seqs
    │   │   └── non_singleton_clusters
    |   |       ├── phylo_singletons.txt
    │   │       └── 1.faa
    │   ├── align
    │   │   └── 1.afa
    │   ├── trim 
    │   │   └── 1.afa
    │   └── phylo
    │       └── 1.treefile
    └── logs

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

### Future development 

#### Containerisation

If you wish to run this in a container (i.e. --use-singularity) please also install [singularity](https://sylabs.io/docs/#singularity)

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity 

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

