# filter_db.py - script to find the correspondence between 
# dbPTM entries and AlphaFold models
# Copyright (C) 2022 Matteo Tiberti, Danish Cancer Society
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import numpy as np
import os


# input
snof = 'S-nitrosylation'
af_models = 'accession_ids.csv'

# output
joined='results/joined_dbs.csv'
filtered='results/filtered_joined_dbs.csv'
selected_entries = "results/selected_proteins.csv"
not_found_fname = "results/not_found_in_af.csv"

if not os.path.exists('./results'):
    os.mkdir('./results')

# reuse time-consuming step if available
# WARNING: This step uses ~80GB of memory!
if os.path.exists(joined):
    print("Will reuse results/filtered_joined_dbs.csv instead of recreating it")

    snos = pd.read_csv(joined)
else:

    snos = pd.read_csv(snof, delim_whitespace=True, header=None, names=['up_id', 'up_ac', 'residue', 'ptm', 'pubmed_id', 'sequence'],
            dtype={'up_id':str, 'up_ac':str, 'residue':str, 'ptm':str, 'pubmed_id':str, 'sequence':str})

    print(f"SNO database contains {len(snos)} entries")
    print(f"These correspond to {snos['up_ac'].unique().shape[0]} unique UniProt ACs")

    af_models = pd.read_csv(af_models, header=None, names=['up_ac', 'first_aa','last_aa', 'fname','version']).set_index('up_ac').sort_index()
    snos = snos.join(af_models, on='up_ac')

    snos.to_csv(joined, index=False)

# filter out proteins for which no match in AF was found
not_found = snos[snos.fname.isna()]
not_found.to_csv(not_found_fname)

snos = snos[~ snos.fname.isna()].reset_index()
print(f"of these, {len(not_found)} were discarded as no corresponding AlphaFold model was found. {len(snos)} are left")

snos.to_csv(filtered, index=False)

snos = pd.read_csv(joined)
snos = snos[~ snos.fname.isna()].reset_index()

# group by Uniprot AC
snos_g = snos.groupby('up_ac').agg(list)

# check if all the entries are defined on the same AF model
snos_g['is_single_structure'] = snos_g.apply(lambda r: len(set(r['fname'])) == 1, axis=1)
assert np.all(snos_g['is_single_structure'])
# Yes! Then we can just consider a 1:1 relationship between Uniprot AC and models

# check that the index doesn't have duplicates
assert snos_g.index.unique().shape == snos_g.index.shape

# collapse multiple same sno to single. Possible as all Uniprot ACs have the same Uniprot IDs
# snos_g['up_id'] = snos_g.apply(lambda r: list(set(r['up_id']))[0], axis=1)
snos.to_csv(selected_entries, index=False)

