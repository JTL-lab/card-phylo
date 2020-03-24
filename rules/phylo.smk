rule align:
    input:
        "seqs/clustered/mmseqs/non_singletons/{clusterid}.faa"
    output:
        "align/{clusterid}.afa"
    log:
        "logs/{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """

rule trim:
    input:
        "align/{clusterid}.afa"
    output:
        "trim/{clusterid}.afa"
    log:
        "logs/{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        trimal -keepheader -in {input} -out {output} -automated1 2> {log}
        """

rule phylo:
    input:
        "trim/{clusterid}.afa"
    output:
        "phylo/{clusterid}.treefile"
    log:
        "phylo/{clusterid}.iqtree"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        iqtree -m LG+G -s {input} -pre phylo/{wildcards.clusterid} -bb 1000
        """


def aggregate_trees(wildcards):
    """
    Aggregate the tree names
    """
    checkpoint_output = checkpoints.write_mmseq_clusters.get(**wildcards).output[0]
    return expand('phylo/{clusterid}.treefile',
                   clusterid=glob_wildcards(os.path.join(checkpoint_output, '{clusterid}.faa')).clusterid)


rule summary:
    input:
        aggregate_trees
    output:
        "pipeline_complete.txt"
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/summarise_phylo.log"
    shell:
        "touch pipeline_complete.txt"
