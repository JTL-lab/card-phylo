rule organise_by_family:
    conda:
        'envs/card-phylo.yml'
    input:
        card_json = f"card/canonical/{config['canonical_version']}/card.json",
        prevalence = f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta"
    output:
        b = "seqs/families/",
        a = "a.txt"
    message:
        "Grouping sequences by annotated family name"
    shell:
        """
        mkdir -p seqs/families
        python scripts/dump_to_gene_family_fasta.py {input.card_json} {input.prevalence} seqs/families
        touch {output.a}
        """

## checkpoint is run for this rule as it has a variable number of outputs
## i.e. number of clusters and we want to re-evaluate 
#rule cluster_by_mmseq:
#    input:
#        "seqs/protein_homogs.fasta"
#    output:
#        "seqs/clustered/mmseqs/amr_clusters_all_seqs.fasta"
#    params:
#        threads = config['threads_per_job'],
#    conda: 
#        "envs/card-phylo.yml"
#    log:
#        "logs/mmseq.log"
#    shell:  
#        """
#        mmseqs easy-cluster --remove-tmp-files 1 --threads {params.threads} {input} seqs/clustered/mmseqs/amr_clusters seqs/clustered/mmseqs/tmp 2>&1 > {log}
#        """
#
#rule write_mmseq_clusters:
#    input:
#        "seqs/clustered/mmseqs/amr_clusters_all_seqs.fasta"
#    output:
#        directory("seqs/clustered/mmseqs/non_singletons")
#    conda: 
#        "../envs/card-phylo.yml"
#    log:
#        "logs/mmseq_cluster_singletons.log"  
#    shell:
#        "python ../scripts/write_mmseqs_clusters.py {input} seqs/clustered/mmseqs/non_singletons > {log}"
#
##rule cluster_by_cdhit:
##    input:
##        "seqs/protein_homogs.fasta"
##    output:
##        "seqs/clustered/cdhit/ asdasdas.fasta"
##    conda:
##        "envs/card-phylo.yml"
##    log:
##        "logs/cdhit.log"
##    shell:
##        """
##        cdhit -i {input} -c 0.5 -o {output} 
##        """
