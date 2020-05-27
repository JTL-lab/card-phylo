rule cdhit_align:
    input:
        "seqs/clustered/cdhit/non_singletons/{clusterid}.faa"
    output:
        "align/cdhit/{clusterid}.afa"
    log:
        "logs/cdhit_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """

rule cdhit_trim:
    input:
        "align/cdhit/{clusterid}.afa"
    output:
        "trim/cdhit/{clusterid}.afa"
    log:
        "logs/cdhit_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        trimal -keepheader -in {input} -out {output} -automated1 2> {log}
        """

rule cdhit_phylo:
    input:
        "trim/cdhit/{clusterid}.afa"
    output:
        "phylo/cdhit/{clusterid}.treefile"
    log:
        "phylo/cdhit/{clusterid}.iqtree"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        iqtree -m LG+G -s {input} -pre phylo/{wildcards.clusterid} -bb 1000
        """


def cdhit_aggregate_trees(wildcards):
    """
    Aggregate the tree names
    """
    checkpoint_output = checkpoints.write_cdhit_clusters.get(**wildcards).output[0]
    return expand('phylo/cdhit/{clusterid}.treefile',
                   clusterid=glob_wildcards(os.path.join(checkpoint_output, '{clusterid}.faa')).clusterid)


rule cdhit_phylo:
    input:
        cdhit_aggregate_trees
    output:
        "cdhit_pipeline_complete.txt"
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/cdhit_phylo.log"
    shell:
        "touch {output}"
