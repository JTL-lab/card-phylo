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
        "pipeline_complete.txt"
    log:
        "logs/pipe.log"

rule mmseqs_phylo_pipe:
    input: 
        "mmseqs_pipeline_complete.txt"

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

rule concatenate_seqs:
    input:
        f"card/canonical/{config['canonical_version']}/protein_fasta_protein_homolog_model.fasta",
        f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta"
    output: 
        "seqs/protein_homogs.fasta"
    log:
        "logs/concatenate.log"
    shell:
        """
        cat {input} > {output}
        sed -i 's/ /_/g' {output}
        echo "Done" > {log}
        """

########################### Cluster Sequences  ##############################

#rule organise_by_family:
#    #conda:
#    #    'envs/card-phylo.yml'
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

rule cluster_by_mmseq:
    input:
        "seqs/protein_homogs.fasta"
    output:
        "seqs/mmseqs_clusters/amr_clusters_all_seqs.fasta"
    params:
        threads = config['threads_per_job'],
    #conda: 
    #    "envs/card-phylo.yml"
    log:
        "logs/mmseq.log"
    shell:  
        """
        mmseqs easy-cluster --remove-tmp-files 1 --threads {params.threads} {input} seqs/mmseqs_clusters/amr_clusters seqs/mmseqs_clusters/tmp 2>&1 > {log}
        """

checkpoint write_mmseqs_clusters:
    input:
        "seqs/mmseqs_clusters/amr_clusters_all_seqs.fasta"
    output:
        cluster_dir = directory("seqs/mmseqs_clusters/non_singleton_clusters")
    log:
        "logs/mmseq_cluster_splitting.log"
    #conda: 
    #    "../envs/card-phylo.yml"
    shell:
        "python ../scripts/write_mmseqs_clusters.py -c {input} -o {output.cluster_dir} 2>&1 > {log}"

rule mmseqs_align:
    input:
        "seqs/mmseqs_clusters/non_singleton_clusters/{clusterid}.faa"
    output:
        "align/mmseqs_{clusterid}.afa"
    log:
        "logs/mmseqs_mafft_{clusterid}.log"
    #conda: 
    #    "envs/card-phylo.yml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """

rule mmseqs_trim:
    input:
        "align/mmseqs_{clusterid}.afa"
    output:
        "trim/mmseqs_{clusterid}.afa"
    log:
        "logs/mmseqs_trimal_{clusterid}.log"
    #conda: 
    #    "envs/card-phylo.yml"
    shell:
        """
        trimal -keepheader -in {input} -out {output} -automated1 2>&1 > {log}
        """

rule mmseqs_tree:
    input:
        "trim/mmseqs_{clusterid}.afa"
    output:
        "phylo/mmseqs_{clusterid}.treefile"
    log:
        "logs/mmseqs_phylo_{clusterid}.log"
    #conda: 
    #    "envs/card-phylo.yml"
    shell:
        """
        iqtree -m LG+G -s {input} -pre phylo/mmseqs_{wildcards.clusterid} -bb 1000 2>&1 > {log}
        """


def aggregate_phylo(wildcards):
    checkpoint_output = checkpoints.write_mmseqs_clusters.get(**wildcards).output[0]
    return expand("phylo/mmseqs_{i}.treefile",
    i=glob_wildcards(os.path.join(checkpoint_output, "{i}.faa")).i)


rule mmseqs_phylo:
    input:
        aggregate_phylo
    output:
        "mmseqs_pipeline_complete.txt"
    #conda: 
    #    "envs/card-phylo.yml"
    shell:
        "touch {output}"
