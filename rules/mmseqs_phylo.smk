rule mmseqs_align:
    input:
        "seqs/clustered/mmseqs/non_singletons/{clusterid}.faa"
    output:
        "align/mmseqs/{clusterid}.afa"
    log:
        "logs/mmseqs_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        mafft --auto {input} > {output} 2> {log}
        """

rule mmseqs_trim:
    input:
        "align/mmseqs/{clusterid}.afa"
    output:
        "trim/mmseqs/{clusterid}.afa"
    log:
        "logs/mmseqs_{clusterid}.log"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        trimal -keepheader -in {input} -out {output} -automated1 2> {log}
        """

rule mmseqs_tree:
    input:
        "trim/mmseqs/{clusterid}.afa"
    output:
        "phylo/mmseqs/{clusterid}.treefile"
    log:
        "phylo/mmseqs/{clusterid}.iqtree"
    conda: 
        "envs/card-phylo.yml"
    shell:
        """
        iqtree -m LG+G -s {input} -pre phylo/{wildcards.clusterid} -bb 1000
        """

rule mmseqs_phylo:
    input:
        dynamic('phylo/mmseqs/{clusterid}.treefile')
    output:
        "mmseqs_pipeline_complete.txt"
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/mmseqs_phylo.log"
    shell:
        "touch {output}"
