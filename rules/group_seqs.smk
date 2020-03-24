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

rule cluster_by_identity:
    input:
        concatenated = "seqs/protein_homogs.fasta"
    output:
        dynamic("seqs/clustered/{clusterid}.fasta")
    params:
        threads = config['threads_per_job'],
    message:
        f"Clustering Protein Sequences with MMSEQS2"
    shell:  
        """
        mkdir -p seqs/clustered/mmseqs/tmp
        mmseqs easy-cluster --threads {params.threads} {input.concatenated} seqs/clustered/mmseqs/amr_clusters seqs/clustered/mmseqs/tmp 2>&1 > /dev/null
        python ../scripts/write_mmseqs_clusters.py seqs/clustered/mmseqs/amr_clusters_all_seqs.fasta seqs/clustered 2>&1 > /dev/null
        """
