#Requirements
- to have generated the data for the mutations of interest with cancermuts in raw_data

#symbolic link to cancermuts metatable
ln -s /data/raw_data/computational_data/cancermuts_data/tcga_3d/p53/pan_cancer_clinvar_others/03102022/metatable_pancancer_TP53.csv

# The script get_mutlist.py takes in input the cancermut's metatable
# and return a mutlist.txt file (mutation format: A556T)
# containing only the mutations with a source and occurring at positions covered by the selected structure

# If you have the full-length protein or a trimmed model in which only the N-ter and/or C-ter portions have been removed
# you can follow the single domain example (i.e., it is enough to specify the N-ter and C-ter residue after the -n and -c flags, respectively).
# If you have trimmed the structure in multiple domains, follow the multiple domains example:
# For instance, if your structures are covering two ranges of residues (e.g., 1:341 and 501:756),
# you have to specify both the N-ter residues space-separated after the -n flag (i.e., -n 1 501) and the C-ter residues after the -c flag (i.e., -c 341 756).

module load python/3.7/modulefile

# single domain
python get_mutlist.py -m metatable_pancancer_TP53.csv -d 91:289 -ch A -M -R -H -C -p afTP53.pdb
#n, c is the residues at the n and c terminal

