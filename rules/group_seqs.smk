#rule organise_by_family:
#    input:
#        card_json = f"card/canonical/{config['canonical_version']}/card.json",
#        prevalence = f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta"
#    output:
#        "seqs/families/dumped.txt"
#    message:
#        "Grouping sequences by annotated family name"
#    shell:
#        """
#        mkdir -p seqs/families
#        python scripts/dump_to_gene_family_fasta.py {input.card_json} {input.prevalence} seqs/families
#        touch seqs/families/dumped.txt
#        """
rule cluster_by_mmseq:
    input:
        "seqs/protein_homogs.fasta"
    output:
        "seqs/clustered/mmseqs/amr_clusters_all_seqs.fasta"
    params:
        threads = config['threads_per_job'],
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/mmseq.log"
    shell:  
        """
        mmseqs easy-cluster --remove-tmp-files 1 --threads {params.threads} {input} seqs/clustered/mmseqs/amr_clusters seqs/clustered/mmseqs/tmp 2>&1 > {log}
        """

checkpoint write_mmseq_clusters:
    input:
        "seqs/clustered/mmseqs/amr_clusters_all_seqs.fasta"
    output:
        directory("seqs/clustered/mmseqs/non_singletons")
    conda: 
        "envs/card-phylo.yml"
    log:
        "logs/mmseq_cluster_singletons.log"  
    shell:
        "python ../scripts/write_mmseqs_clusters.py {input} seqs/clustered/mmseqs/non_singletons > {log}"
