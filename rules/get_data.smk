# An example collection of Snakemake rules imported in the main Snakefile.
#

rule download_card:
    output:
    
    script:
        """
        scripts/get_data.sh
        """

rule get_card:
    output: 
        
        
        
   

rule organise_by_family:


rule cluster_by_identity:
    
      
