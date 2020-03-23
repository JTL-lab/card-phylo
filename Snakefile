# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
configfile: "config.yaml"
report: "report/workflow.rst"

# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

rule all:
    input:
        "phlyo/{amr}.tree"

rule align:
    input:
        "seqs/{amr}.fna"
    output:
        "align/{amr}.afa"
    shell:
        "kalign < {input} > {output} 2> /dev/null"

rule trim:
    input:
        "align/{amr}.afa"
    output:
        "trim/{amr}.afa"
    shell:
        "trimal -in {input} -out {output} -automated1 2> /dev/null"

rule phylo:
    input:
        "trim/{amr}.afa"
    output:
        "phylo/{amr}.tree"
    shell:
        "fasttree -gtr -nt {input} > {output} 2> /dev/null"

include: "rules/get_card.smk"
