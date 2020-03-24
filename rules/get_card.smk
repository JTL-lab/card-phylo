# An example collection of Snakemake rules imported in the main Snakefile.
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
    conda: 
        "envs/card-phylo.yml"
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
    conda:
        "envs/card-phylo.yml"
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
    conda:
        "envs/card-phylo.yml"
    log:
        "logs/concatenate.log"
    shell:
        """
        cat {input} > {output}
        sed -i 's/ /_/g' {output}
        echo "Done" > {log}
        """
