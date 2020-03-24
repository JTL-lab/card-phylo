# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
from glob import glob

configfile: "config.yaml"
report: "report/workflow.rst"
workdir: f"card-phylo_canonical_{config['canonical_version']}_prevalence_{config['prevalence_version']}"

rule all:
    input:
        "done"
#report/report.rst"

rule report:
    input:
        "phlyo/{amr}.tree"
    output:
        "report/report.rst"

include: "rules/get_card.smk"
include: "rules/group_seqs.smk"
include: "rules/phylo.smk"
