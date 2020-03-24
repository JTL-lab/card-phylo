rule align:
    input:
        dynamic("seqs/clustered/{clusterid}.fasta")
    output:
        "align/{clusterid}.afa"
    log:
        "logs/{clusterid}.log"
    message:
        "Aligning sequences"
    shell:
        """
        mkdir -p align
        mafft --auto {input} > {output} 2>&1 >> {log}
        """

rule trim:
    input:
        "align/{clusterid}.afa"
    output:
        "trim/{clusterid}.afa"
    log:
        "logs/{clusterid}.log"
    message:
        "Trimming alignments"
    shell:
        """
        mkdir -p trim
        trimal -in {input} -out {output} -automated1 2>&1 >> {log}"
        """

rule phylo:
    input:
        "trim/{clusterid}.afa"
    output:
        "phylo/{clusterid}.tree"
    log:
        "logs/{clusterid}.log"
    message:
        "Inferring ML Phylogeny"
    shell:
        """
        mkdir -p phylo
        fasttree {input} > {output} 2>&1 >> {log}
        """

rule summary:
    input:
        dynamic("phylo/{clusterid}.tree")
    output:
        "done"
    shell:
        "touch done"
