# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
configfile: "config.yaml"
report: "report/workflow.rst"
workdir: f"card-phylo_canonical_{config['canonical_version']}_prevalence_{config['prevalence_version']}"

#rule all:
#    input:
#        "done"
#report/report.rst"

rule all:
    input:
        "mmseqs_complete.txt"

############################## Download and Parse CARD Databases ##############################
rule extract_card:
    input:
        canonical = f"card/canonical/{config['canonical_version']}/broadstreet-v{config['canonical_version']}.tar.bz2",
        prevalence = f"card/prevalence/{config['prevalence_version']}/prevalence-v{config['prevalence_version']}.tar.bz2"
    output:
        f"card/canonical/{config['canonical_version']}/card.json",
        f"card/canonical/{config['canonical_version']}/protein_fasta_protein_homolog_model.fasta",
        f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta",
    params:
        canonical_version = config["canonical_version"],
        prevalence_version = config["prevalence_version"]
    log:
        "logs/card_extract.log"
    shell:
        """
        tar -C card/canonical/{params.canonical_version} -xvf {input.canonical} 2>&1 >> {log}
        tar -C card/prevalence/{params.prevalence_version} -xvf {input.prevalence} 2>&1 >> {log}
        gunzip card/prevalence/{params.prevalence_version}/protein_fasta_protein_homolog_model_variants.fasta
        """

rule download_card:
    output:
        f"card/canonical/{config['canonical_version']}/broadstreet-v{config['canonical_version']}.tar.bz2",
        f"card/prevalence/{config['prevalence_version']}/prevalence-v{config['prevalence_version']}.tar.bz2",
    params:
        canonical_version = config["canonical_version"],
        prevalence_version = config["prevalence_version"]
    log:
        "logs/card_download.log"
    shell:
        """
        wget -P card/canonical/{params.canonical_version} https://card.mcmaster.ca/download/0/broadstreet-v{params.canonical_version}.tar.bz2 2>&1 >> {log}
        wget -P card/prevalence/{params.prevalence_version} https://card.mcmaster.ca/download/6/prevalence-v{params.prevalence_version}.tar.bz2 2>&1 >> {log}
        """

# only include canonical for testing purposes
rule concatenate_seqs:
    input:
        f"card/canonical/{config['canonical_version']}/protein_fasta_protein_homolog_model.fasta",
        #f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta"
    output: 
        "card_protein_homogs.fasta"
    log:
        "logs/concatenate.log"
    shell:
        """
        cat {input} > {output}
        sed -i 's/ /_/g' {output}
        echo "Done" > {log}
        """

########################### Phylo Generation ###################################

rule align:
    input:
        "{clustertype}/seqs/non_singleton_clusters/{clusterid}.faa"
    output:
        "{clustertype}/align/{clusterid}.afa"
    log:
        "logs/{clustertype}/mafft_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """

rule trim:
    input:
        "{clustertype}/align/{clusterid}.afa"
    output:
        "{clustertype}/trim/{clusterid}.afa"
    log:
        "logs/{clustertype}/trimal_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        trimal -keepheader -in {input} -out {output} -automated1 2>&1 > {log}
        """

rule tree:
    input:
        "{clustertype}/trim/{clusterid}.afa"
    output:
        "{clustertype}/phylo/{clusterid}.treefile"
    log:
        "logs/{clustertype}/phylo_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    params:
        prefix = "{clustertype}/phylo/{clusterid}"
    shell:
        """
        iqtree -fast -m LG+G -s {input} -pre {params.prefix} 2>&1 > {log}
        """

######################### MMSeqs Clustering ####################################

rule cluster_by_mmseq:
    input:
        "card_protein_homogs.fasta"
    output:
        "mmseqs/seqs/amr_clusters_all_seqs.fasta"
    params:
        threads = config['threads_per_job'],
        prefix = "mmseqs/seqs/amr_clusters",
        tmp = "mmseqs/seqs/mmseqs_clusters/tmp"
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/mmseqs/mmseq_cluster.log"
    shell:  
        """
        mkdir -p {params.tmp}
        mmseqs easy-cluster --remove-tmp-files 1 --threads {params.threads} {input} {params.prefix} {params.tmp} 2>&1 > {log}
        """

checkpoint write_mmseqs_clusters:
    input:
        "mmseqs/seqs/amr_clusters_all_seqs.fasta"
    output:
        cluster_dir = directory("mmseqs/seqs/non_singleton_clusters")
    log:
        "logs/mmseqs/mmseq_cluster_splitting.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        "python ../scripts/write_mmseqs_clusters.py -c {input} -o {output.cluster_dir} 2>&1 > {log}"

def aggregate_mmseqs_phylo(wildcards):
    checkpoint_output = checkpoints.write_mmseqs_clusters.get(**wildcards).output[0]
    return expand("mmseqs/phylo/{i}.treefile",
    i=glob_wildcards(os.path.join(checkpoint_output, "{i}.faa")).i)

rule mmseqs:
    input:
        aggregate_mmseqs_phylo
    output:
        "mmseqs_complete.txt"
    shell:
        "touch {output}"


###################### Family Grouping ########################################

#rule family:
#    output

#rule mmseqs:
#    input: 

#rule id70:
#
#rule id75:
#
#rule id80:
#
#rule id85:
#
#rule id90:
#
#rule id95:

#rule organise_by_family:
#    conda:
#        'envs/card-phylo.yml'
#    input:
#        card_json = f"card/canonical/{config['canonical_version']}/card.json",
#        prevalence = f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta"
#    output:
#        "seqs/families/",
#    message:
#        "Grouping sequences by annotated family name"
#    shell:
#        """
#        mkdir -p seqs/families
#        python ../scripts/dump_to_gene_family_fasta.py -c {input.card_json} -p {input.prevalence} -o {output} 
#        """

#include "rules/mmseqs_phylo.smk"
#include "rules/family_phylo.smk"
#include "
#
