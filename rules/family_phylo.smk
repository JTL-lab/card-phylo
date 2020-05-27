rule family_align:
    input:
        "seqs/family/non_singletons/{clusterid}.faa"
    output:
        "align/family/{clusterid}.afa"
    log:
        "logs/family_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """

rule family_trim:
    input:
        "align/family/{clusterid}.afa"
    output:
        "trim/family/{clusterid}.afa"
    log:
        "logs/family_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        trimal -keepheader -in {input} -out {output} -automated1 2> {log}
        """

rule family_phylo:
    input:
        "trim/family/{clusterid}.afa"
    output:
        "phylo/family/{clusterid}.treefile"
    log:
        "phylo/family/{clusterid}.iqtree"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        iqtree -m LG+G -s {input} -pre phylo/{wildcards.clusterid} -bb 1000
        """


def family_aggregate_trees(wildcards):
    """
    Aggregate the tree names
    """
    checkpoint_output = checkpoints.write_mmseq_clusters.get(**wildcards).output[0]
    return expand('phylo/family/{clusterid}.treefile',
                   clusterid=glob_wildcards(os.path.join(checkpoint_output, '{clusterid}.faa')).clusterid)


rule family_phylo:
    input:
        family_aggregate_trees
    output:
        "family_pipeline_complete.txt"
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/family_phylo.log"
    shell:
        "touch {output}"
