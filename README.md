# UMRAD
Universal Multi-omics Reference and Alignment Database

## What is it?
This repository involves a series of pipelines that create a Universal Reference to unite and analyze multi-omics data.
Currently the different databases are separated - to keep things simple for me when I made them.
    <br>1. Universal Taxonomy Database: [found here](https://github.com/TealFurnholm/Universal-Taxonomy-Database)
    <br>2. Universal Compounds Database and + 3. Universal Reactions Database: [found here](https://github.com/TealFurnholm/Universal_Biological_Compounds_Database)
    <br>4. Universal Protein Alignment Database: [found here](https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database)
    <br>5. Universal ncRNA Alignment Database: [found here](https://github.com/TealFurnholm/Fix_RNACentral_Taxonomy)
<p>

## Why Universal?
- These databases span all kingdoms of life. 
- The databases allow the simultaneous identification of microbial community phylogeny and functions. 
- All the biological molecules/compounds and their data are linked to their enzymes and transporters to map the flow of metabolites in microbe-microbe or microbe-host interactions. 
- The databases are used for (meta)transcriptomics, (meta)proteomics, metabolomics, metagenomics, and for novel binning and MAG quality control software I've created. 
This way many types of data can be combined into a secondary analysis, with taxonomy and functions being directly linked. It also covers both the protein and the non-coding fraction of sequencing. 

## How to Use
There are several primary and secondary microbiome analysis pipelines on my main GitHub page that use these databases:
* RNAseq/Metatranscriptome Analysis [here](https://github.com/TealFurnholm/Strain-Level_Metatranscriptome_Analysis)
* Metagenome Primary Analysis [here](https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis)
* Metabolomics Analysis - TBD; for now you can directly link the metabolomics output to compounds->reactions->proteins/oranisms using the Functional and Protein databases
* Strain-level Metagenome Binning [here](https://github.com/TealFurnholm/Community-Based-Metagenome-Binning)
* Secondary Functional Analysis and Visualization [here](https://github.com/TealFurnholm/Meta-omics_Functional_Analysis)    
    
## Start: [here](https://github.com/TealFurnholm/Universal-Taxonomy-Database)
