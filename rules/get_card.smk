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
    message:
        "Extracting CARD databases"
    log:
        "run.log"
    shell:
        """
        tar -C card/canonical/{params.canonical_version} -xvf {input.canonical} 2>&1 >>  {log}
        tar -C card/prevalence/{params.prevalence_version} -xvf {input.prevalence} 2>&1 >>  {log}
        gunzip card/prevalence/{params.prevalence_version}/protein_fasta_protein_homolog_model_variants.fasta
        """

rule download_card:
    output:
        f"card/canonical/{config['canonical_version']}/broadstreet-v{config['canonical_version']}.tar.bz2",
        f"card/prevalence/{config['prevalence_version']}/prevalence-v{config['prevalence_version']}.tar.bz2",
    log: 
        "run.log"
    message:
        f"Downloading CARD Canonical ({config['canonical_version']}) and Prevalence datasets ({config['prevalence_version']})"
    params:
        canonical_version = config["canonical_version"],
        prevalence_version = config["prevalence_version"]
    shell:
        """
        mkdir -p card/{{canonical,prevalence}}
        mkdir -p card/canonical/{params.canonical_version}
        wget -P card/canonical/{params.canonical_version} https://card.mcmaster.ca/download/0/broadstreet-v{params.canonical_version}.tar.bz2 2> {log}
        mkdir -p card/prevalence/{params.prevalence_version}
        wget -P card/prevalence/{params.prevalence_version} https://card.mcmaster.ca/download/6/prevalence-v{params.prevalence_version}.tar.bz2 2>&1 >> {log}
        """

rule concatenate_seqs:
    input:
        canonical  = f"card/canonical/{config['canonical_version']}/protein_fasta_protein_homolog_model.fasta",
        prevalence = f"card/prevalence/{config['prevalence_version']}/protein_fasta_protein_homolog_model_variants.fasta"
    output: 
        concatenated = "seqs/protein_homogs.fasta"
    message:
        "Concatenating all CARD protein fasta homolog sequences and removing spaces"
    shell:
        """
        mkdir -p "seqs"
        cat {input.canonical} {input.prevalence} > {output.concatenated}
        sed -i 's/ /_/g' {output.concatenated}
        """


