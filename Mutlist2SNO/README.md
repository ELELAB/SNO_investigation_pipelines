# Introduction

Mutlist2SNO pipeline is a workflow to filter list of mutations based on the integration of information collected in the input file for the SNO_model pipeline and output files from the SNO_finder pipeline. In particular, it takes as input:
-pdb files of SNO containing proteins, 
-the mavisp csv files for each protein of interest (containing the information about the mutations retrieved from several databases in terms of stability),
-the input csv file of SNO_finder pipeline containig information about the proteins of interest.
-the output file from SNO_finder pipeline containing the proximal and the SNO cysteines for each proteins of interest. 
The pipeline runs PyInteraph on the pdb file selecting those residues in the around of 8 A° establishing contacts with a proximal or a SNO cysteine. Subsequently, it returns as output two files containing the Neutral and Uncertain mutations in terms of ddg of stability on the selected residues expressed with one letter annotation or in the cabsflex format.

The pipeline:

 1. Runs the PyInteraph2 software on the pdb file provided as input in the appropiate folder for each protein to calculate cmPSN with 8A distance cut-off from a proximal cysteine or a cysteine forming a SNO site, which is meant to identify contacts between each pair of residues in the protein
 (results/pynteraph/{gene_name}/{uniprot_ac}_{start}-{end}_cmpsn_all.csv).

 2. Examines the input file for SNO_model pipeline and the output file from SNO_finder pipeline (proteins.csv and cys_database_far_seq_acc10_human.csv respectively) in order to detect the proximal and SNO cysteines for each protein (starting_dataset.csv).

 3. Combines the output obtained in the step 1 and the information on the cysteines from the step 2, in order to retrieve only the residues establishing contacts with the cysteine of interest (results/pynteraph/{gene_name}/{uniprot_ac}_{start}-{end}_cmpsn_all_parsed.csv).

 4. Compares the mavisp.csv file provided as input ({GENE_NAME}-simple_mode.csv) with the residues collected in the output obtaind in the step 3 in order to select those mutation with a Neutral or an Uncertain effect in terms of protein stability, affecting those reesidues selected in the step3.

 5. Annotates the mutations retrieved in step 4 in two ways:  
    - creates a mutlist.txt (subset_{uniprot_ac}_{start}-{end}_mutlist.txt) file containing the mutations expressed with one letter code ([A-Z][0-9][A-Z]) 
    - creates a cabsflex_mutlist.txt (subset_{uniprot_ac}_{start}-{end}_mutlist_cabsflex.txt) containing the mutations expressed accordingly to the cabsflex format (pdbfile, chain, res_position, res_WT, res_mut)

Since the pipeline exploits PyInteraph to compute the contacts, the Glycines will be not considered. So, the mutations regarding the Glycines in the around of 8 A° from a proximal or a SNO cysteine will be not retrieved.

# Requirements

In order to replicate the analysis you will need to have the following software.
Versions the pipeline has been tested with and used are indicated in 
parenthesis.

  - Python >=3.7 with packages 
      re
      pandas (1.3.2)
      os
      snakemake (7.17.1)
  - PyInteraph2 (1.1)

# Input and output structures

Input:

The input structure should be like the one in the example below.
The essential files are:
-mutlist_cabsflex.txt, 
-mutlist.txt, 
-csv file form mavisp database (for example CALR-simple_mode.csv)
-pdb_file of the protein of interest (for example P27797_17-404.pdb, the pdb file must be named as follow [uniprot_ac_id]_[start_residues]-[end_residue].pdb).

The other files are not essential but they will be present if the folders (cancermuts, structure selection) for each protein are gotten from the mavisp workflow outputs.


calr
├── cancermuts
│   ├── mutlist_cabsflex.txt
│   ├── mutlist.txt
│   └── readme.txt
├── mavisp_csv
│   └── CALR-simple_mode.csv
└── structure_selection
    └── alphafold_db
        ├── calr
        │   ├── P27797_17-404.pdb
        │   ├── P27797.csv
        │   ├── P27797_dssp.out
        │   ├── P27797.json
        │   └── P27797.pdb
        ├── config_alphafolddb.yaml
        ├── get_alphafolddb_data.py
        ├── readme.txt
        └── regions_pLDDT.csv

Output:

results
├── mutlists
│   ├── calr
│       ├── subset_P27797_17-404_mutlist_cabsflex.txt
│       └── subset_P27797_17-404_mutlist.txt
└── pynteraph
    ├── calr
        ├── P27797_17-404_cmpsn_all.csv
        └── P27797_17-404_cmpsn_all_parsed.csv

# Performing the analysis:

1. Activate the appropriate Python environment if needed

module load python/3.7/modulefile

2. run the pipeline:

snakemake -s Snakefile -c 1


