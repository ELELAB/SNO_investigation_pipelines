#Requirements
- to have generated the data for the mutations of interest with cancermuts in raw_data

#symbolic link to cancermuts metatable
ln -s /data/raw_data/computational_data/cancermuts_data/sno/GRIN2A/pancancer_clinvar/16112022/metatable_pancancer_GRIN2A.csv

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
python get_mutlists.py -m metatable_pancancer_GRIN2A.csv -d 32:387 -ch A -M -R

# multiple domains 
python get_mutlist.py -m metatable_pancancer_MLH1.csv -n 1 501 -c 341 756 
