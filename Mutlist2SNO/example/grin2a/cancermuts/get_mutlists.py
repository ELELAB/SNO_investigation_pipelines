#!/usr/bin/env python
"""
@author: Ludovica Beltrame
"""
import pandas as pd
import argparse
from Bio.SeqUtils import seq3

# Add arguments required to the script to a argparse.ArgumentParser instance.
description = "Get different mutation lists from cancermuts metatable"
parser = argparse.ArgumentParser(description = description)

m_helpstr = "Input metatable"
parser.add_argument("-m", "--metatable", action="store", type=str, help=m_helpstr, required=True)

d_helpstr = "Range of each domain"
parser.add_argument("-d", "--domain", nargs='+', type=str, help=d_helpstr, required=True)

ch_helpstr = "chain of interest"
parser.add_argument("-ch", "--chain", type=str, help=ch_helpstr)

p_helpstr = "pdb file of interest"
parser.add_argument("-p", "--pdbfile", type=str, help=p_helpstr)

mutatex_helpstr = "generate mutatex mutlist"
parser.add_argument("-M", "--mutatex", action='store_true', help=mutatex_helpstr)

rosetta_helpstr = "generate rosetta mutlist"
parser.add_argument("-R", "--rosetta", action='store_true', help=rosetta_helpstr)

hgvs_helpstr = "generate HGVS mutlist"
parser.add_argument("-H", "--hgvs", action='store_true', help=hgvs_helpstr)

cabsflex_helpstr = "generate cabsflex mutlist"
parser.add_argument("-C", "--cabsflex", action='store_true', help=cabsflex_helpstr)

args = parser.parse_args()

def oneletter_format(dataframe):
    '''Mutation format: wt position mut (e.g., A75C)'''
    
    #Add mutation column
    dataframe['mutation'] = dataframe['ref_aa'] + dataframe["aa_position"].astype(str) + dataframe["alt_aa"]
    
    return dataframe['mutation'].to_csv('mutlist.txt', header=None, index=None)

def mutatex_format(dataframe, chain):
    '''Mutation format: wt chain position mut (e.g., AA75C)'''

    #Add mutation column
    dataframe['mutatex_mutation'] = dataframe['ref_aa'] + chain + dataframe["aa_position"].astype(str) + dataframe["alt_aa"]
    
    return dataframe['mutatex_mutation'].to_csv('mutlist_mutatex.txt', header=None, index=None)

def mutatexP_format(dataframe, chain):
    '''Mutation format: wt chain position mut (e.g. TA75C)
    and wt chain position phospho_type (e.g., TA75p)'''

    #Define dictionary of phosphorylations
    phospho_dic = {'T': 'p', 'S': 's', 'Y': 'y'}
    
    #Filter dataframe for positions supporting phosphorylation
    phospho_dataframe = dataframe[dataframe['phosphorylation_site'] == 'P'].copy()
    
    #Add mutations column
    phospho_dataframe['mutatex_mutation'] = phospho_dataframe['ref_aa'] + chain + phospho_dataframe["aa_position"].astype(str) + phospho_dataframe["alt_aa"]
    phospho_dataframe['phospho_type'] = phospho_dataframe['ref_aa'].apply(lambda x: phospho_dic.get(x))
    phospho_dataframe['mutatex_P_mutation'] = phospho_dataframe['ref_aa'] + chain + phospho_dataframe["aa_position"].astype(str) + phospho_dataframe['phospho_type']
    
    return pd.concat([phospho_dataframe["mutatex_mutation"],phospho_dataframe["mutatex_P_mutation"].drop_duplicates()], axis = 0).to_csv('mutlist_mutatex_P.txt', header=None, index=None)

def rosetta_format(dataframe, chain):
    '''Mutation format: chain.wt.position.mut  chain (e.g., A.A.75.C A)'''

    #Add mutation column
    dataframe['rosetta_mutation'] = chain + '.' + dataframe['ref_aa'] + '.' + dataframe["aa_position"].astype(str) + '.' + dataframe["alt_aa"] + ' ' + chain
    
    return dataframe['rosetta_mutation'].to_csv('mutlist_rosetta.txt', header=None, index=None)

def hgvs_format(dataframe):
    '''Mutation format: p. wt3L position mut3L (e.g., p.Ala75Cys)'''

    #Add mutation column
    dataframe['mutation_hgvs'] = 'p.' + dataframe['wt_3L'] + dataframe["aa_position"].astype(str) + dataframe["mut_3L"]

    return dataframe['mutation_hgvs'].to_csv('mutlist_hgvs.txt', header=None, index=None)

def cabsflex_format(dataframe, pdb, chain):
    '''Mutation format: pdb chain position wt mut (e.g., 2XWR.pdb A 75 ALA CYS)'''

    #Add mutation column
    dataframe['mutation_cabsflex'] = pdb + ' ' + chain + ' ' + dataframe['aa_position'].astype(str) + ' ' + dataframe['wt_3L'].str.upper() + ' ' + dataframe['mut_3L'].str.upper()

    return dataframe['mutation_cabsflex'].to_csv('mutlist_cabsflex.txt', header=None, index=None)

if __name__ == '__main__':
    
    cancermuts_metatable = args.metatable
    domain_range = args.domain
    chainID = args.chain
    pdbID = args.pdbfile
    
    # Read cancermuts metatable
    df = pd.read_csv(cancermuts_metatable, index_col=0)
    df.shape

    # Filter out lines without source
    df = df[pd.notnull(df['sources'])]
    df.shape

    # Filter out lines with aa_position not covered by the structure
    filtered_df = pd.DataFrame()
    for r in domain_range:
        n, c = r.split(':')
        print('Nter: ', n)
        print('Cter: ', c)
        filtered_df = filtered_df.append(df[df.aa_position.between(left=int(n), right=int(c), inclusive='both')])
    
    #Convert letter code of wt and mutant residues
    filtered_df['wt_3L'] = filtered_df['ref_aa'].apply(seq3)
    filtered_df['mut_3L'] = filtered_df['alt_aa'].apply(seq3) 
    
    #Generate one letter mutlist 
    oneletter_format(filtered_df)
    
    #Generate mutatex mutlist
    if args.mutatex:
        if chainID:
            mutatex_format(filtered_df, chainID)
            mutatexP_format(filtered_df, chainID)
        else:
            print('Missing chain ID')

    #Generate rosetta mutlist
    if args.rosetta:
        if chainID:
            rosetta_format(filtered_df, chainID)
        else:
            print('Missing chain ID')
    
    #Generate HGV mutlist
    if args.hgvs:
        hgvs_format(filtered_df)
    
    #Generate cabsflex mutlist
    if args.cabsflex:
        if not chainID:
            print('Missing chain ID')
        elif not pdbID:
            print('Missing PDB ID')
        else:
            cabsflex_format(filtered_df, pdbID, chainID)

