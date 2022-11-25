# Introduction

this Snakemake pipeline starts from a list of curate S-nytrosilation sites in several
proteins from different organisms and tries to identify Cysteine residues close to potential
S-Nytrosilation sites on AlphaFold2 models of them. Two residues are considered close if they
are in contact in the cmPSN analysis of PyInteraph2, using 8 A as a distance cut-off
(https://github.com/ELELAB/pyinteraph2/). It also annotates both the SNO site the
found cysteine residues with their pLDDT score in the model and relative solvent
accessible surface (%) of their side-chain as calculated by naccess (http://www.bioinf.manchester.ac.uk/naccess/)

The list of S-nytrosilation sites was obtained from dbPTM (https://awi.cuhk.edu.cn/dbPTM/), Linux/Mac version.
The file was then manually modified to fix header and footer which contained spurious binary data.

This pipeline also needs access to the AlphaFold Protein Structure Database list of models,
in order to be able to tell which protein's models are available.
(https://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv)

The pipeline:
  1. Filters the SNO database file to keep only proteins whose models are available.
  This is done by running the filter_db.py script. Be careful - this script uses up to
  80GB of memory! It writes a results/joined_dbs.csv which adds information
  on the AF2 available models on the original dataset, a results/selected_proteins.csv
  which contains residues and proteins for which a model was available and
  not_found_in_af.csv which contains residues and proteins for which a model was
  not available
  2. For each protein and SNO site, it performs the following operations in
  a separate folder (results/{uniprot_id}):
    1. Downloads the corresponding AlphaFold model (results/{uniprot_id}/{uniprot_ac}.pdb)
    2. Writes a file with corresponding per-residue pLDDT scores (results/{uniprot_id}/{uniprot_ac}.csv)
    3. Runs the PyInteraph software to calculate cmPSN with 8A distance cut-off, which
    is meant to identify contacts between each pair of residues in the protein
    (results/{uniprot_id}/{uniprot_ac}_cmpsn_all.csv)
    4. Runs naccess on the model PDB file. Please notice that we used a modified
    version of naccess in which the maximum number of cubes was increased to allow
    support for larger proteins. (results/{uniprot_id}/{uniprot_ac}.rsa)
    5. From the PyInteraph results, it identifies all the Cysteine residue that are
    in contact with the SNO-Cys in the database. It also collects
       * pLDDT for SNO site and proximal cysteine in the AlphaFold model
       * secondary structure definition for SNO site and proximal cysteine
         according to DSSP
       * solvent accessible surface of side-chain of both the SNO site and
         proximal cysteine, calculated with NACCESS
       * distance (in A) between the SG atoms of the side-chains of SNO site
         and proximal cysteine
       * predicted pKa for both the SNO site and proximal cysteine, calculated
         with PropKa
    These results are written in a single dataset (results/{uniprot_id}/{uniprot_ac}_potential_sno.csv)
  3. aggregates back the information on all the found sites in each dataset and filters out cases for
     which the SNO site identified in dbPTM was actually not a cysteine residue in the AlphaFold model
     (these are saved in results/cys_database_no_SNO_cys.csv). The resulting dataset is saved as
     results/cys_database.csv. Of the identified SNO sites/cysteine pairs, it also saves separately
     those cases for which the sequence distance between SNO site and proximal cysteine is <= 8
     (results/cys_database_close_seq) or > 8 A (results/cys_database_far_seq.csv). Finally, it filters
     the latter dataset again by keeping those cases for which the protein is human and the
     solvent accessible surface of the SNO site residue is >= 10%

# Requirements

In order to replicate the analysis you will need to have the following software.
Versions the pipeline has been tested with and used are indicated in 
parenthesis.

  - Python >=3.7 with packages 
      numpy (1.21.2)
      pandas (1.3.2)
      biopython (1.79)
      snakemake (7.17.1)
      MDAnalysis (1.0.0)
  - PyInteraph2 (1.1)
  - naccess (2.1.1)
  - dssp (3.0.10)
  - propKa (3.4.0)

# Performing the analysis

1. Be sure to have linked or downloaded the AlphaFold structure database
list of models (https://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv)
and the S-nytrosilation database from PTMdb in the current directory.
They should be called accession_ids.csv and S-nytrosilation respectively.

2. Activate the appropriate Python environment if needed

. /usr/local/envs/py37/bin/activate

3. Run the Python script to perform step 1 of the pipeline. 
WARNING: this will use a large amount of memory (up to 80GB)

python filter_db.py

4. run all the other steps of the pipeline:

snakemake -c 12
