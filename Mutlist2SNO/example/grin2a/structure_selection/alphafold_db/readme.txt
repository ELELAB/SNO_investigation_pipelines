# Activate the appropriate Python environment

. /usr/local/envs/py37/bin/activate

# Run the Python script to:
# - Download the AlphaFold structures from the AlphaFold
#   Protein Structure Database.
# - Get the protein regions associated with high/low pLDDT scores.
# - Run DSSP to get the secondary structure assigned to each
#   residue.

python get_alphafolddb_data.py -c config_alphafolddb.yaml
