#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

# get_alphafolddb_data.py - script to download AlphaFold models and their
# associated metadata
# Copyright (C) 2022 Valentina Sora, Danish Cancer Society
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


# Standard library
import itertools
import os
import subprocess
# Third-party packages
import Bio.PDB as bpdb
import pandas as pd
import requests
import yaml



class DSSPParser:

    """Class implementing a parser for DSSP output files.

    This class is taken from the `dssp_tools.py` module
    of the `SLiMfast` Python package.

    `SLiMfast` can be found at https://github.com/ELELAB/SLiMfast.
    """

    # Standard mapping of the DSSP symbols to the symbols
    # that will be used in the FASTA-like DSSP sequence
    # (default mapping is each symbol to itself)
    DSSP_MAPPING = \
        {"G" : "G", "I" : "I", "H" : "H", "B" : "B",
         "S" : "S", "E" : "E", "T" : "T"}

    # Rrepsentation of white spaces (= disordered residues)
    # in the FASTA-like DSSP sequence
    WHITESPACE_REPR = "-"

    # Data types for each column of the final dataframe
    DTYPES = {\
        "RESIDUE" : "int64",
        "CHAIN" : "object",
        "AA" : "object",
        "STRUCTURE" : "object",
        "3-TURNS_HELIX" : "object",
        "4-TURNS_HELIX" : "object",
        "5-TURNS_HELIX" : "object",
        "GEOMETRICAL_BEND" : "object",
        "CHIRALITY" : "object",
        "BETA_BRIDGE_LABEL" : "object",
        "BP1" : "int64",
        "BP2" : "object",
        "ACC" : "int64",
        "N-H-->O_1" : "object",
        "O-->H-N_1" : "object",
        "N-H-->O_2" : "object",
        "O-->H-N_2" : "object",
        "TCO" : "float64",
        "KAPPA" : "float64",
        "ALPHA" : "float64",
        "PHI" : "float64",
        "PSI" : "float64",
        "X-CA" : "float64",
        "Y-CA" : "float64",
        "Z-CA" : "float64"}


    def __init__(self,
                 dssp_mapping = None,
                 whitespace_repr = None):
        """Initialize the DSSP output parser.
        
        Parameters:
        -----------
        dssp_mapping : `dict` or `None`, default: `None`
            Dictionary mapping the 7 symbols of the DSSP dictionary
            to custom symbols.
        whitespace_repr : `str`, default: `"-"`
            How the whitespaces (assigned by DSSP to unstructured
            residues) are represented in the FASTA-like representation
            of the DSSP sequence.
        """

        # Get the mapping for DSSP symbols
        dssp_mapping = \
            dssp_mapping if dssp_mapping is not None \
            else self.DSSP_MAPPING

        # Get the whitespace representation
        whitespace_repr =  \
            whitespace_repr if whitespace_repr is not None \
            else self.WHITESPACE_REPR

        # Add the mapping for whitespaces to the final mapping that
        # will be used
        self.mapping = {**dssp_mapping, " " : whitespace_repr}


    def parse(self,
              dssp_output,
              verbose = False):
        """Parse the output of DSSP and return a dataframe
        with the data.
        Parameters
        ----------
        dssp_output : `str`
            File where the DSSP output is stored.
        verbose : `bool`, default: `False`
            Whether the DSSP output has been generated in verbose
            mode.
        """

        with open(dssp_output, "r") as f:       

            # If 'verbose' has been set to True
            if verbose:
                # Raise an exception, since parsing for output files
                # generated in verbose mode has not been implemented
                # yet
                errstr = \
                    "Only parsing for non-verbose outputs has been " \
                    "implemented so far."
                raise NotImplementedError(errstr)

            # Create an empty list to store data from the DSSP file
            # that will be later converted into a dataframe
            data = []

            # Create an empty dictionary to store the FASTA sequence
            # (per chain) of the structure on which DSSP has been run
            seqdssp = {}

            # Create an empty dictionary to store the FASTA-like DSSP
            # sequence (per chain) of the structure on which DSSP
            # has been run
            fasta = {}

            # Initialize the boolean switch that indicates when
            # to start parsing to False
            start_parse = False

            # Initialize the string that will temporarily hold the
            # FASTA sequence for a single chain to an empty string
            tmp_fasta = ""

            # Initialize the string that will temporarily hold the
            # FASTA-like DSSP sequence for a single chain to an
            # empty string 
            tmp_seqdssp = ""

            # Initialize the string that will temporarily hold the
            # ID of the chain to which the residue prior to the one
            # currently parsed belong to an empty string
            prev_chain = ""

            # Initialize the string that will temporarily hold the
            # residue prior to the one currently parsed to an
            # empty string
            prev_res = ""

            for line in f:

                # Get rid of the newline delimiter
                line = line.rstrip("\n")

                # Start of the DSSP data
                if line[2] == "#":
                    # Turn on the boolean switch for parsing
                    start_parse = True
                    # Skip the current line
                    continue

                # If the boolean switch is off
                if not start_parse:
                    # Skip the current line
                    continue

                # Get residue number, chain ID and residue type
                residue = line[6:10].lstrip(" ")
                chain = line[11].lstrip(" ")
                aa = line[13:15].rstrip(" ")

                # If you reached a chain break
                if aa == "!*":
                    # Save the FASTA-like DSSP sequence of
                    # the chain that just ended
                    seqdssp[prev_chain] = tmp_seqdssp
                    # Save the FASTA sequence of the chain
                    # that just ended
                    fasta[prev_chain] = tmp_fasta
                    # Re-initialize the temporary holders of
                    # the FASTA sequence and of the FASTA-like
                    # DSSP sequence
                    tmp_fasta = ""
                    tmp_seqdssp = ""
                    # Go to the next line
                    continue
                
                # If you found a missing residue
                if aa == "!":
                    # Warn the user
                    warnstr = \
                        f"Missing residue at position " \
                        f"{int(prev_res)+1} chain {prev_chain}."
                    logger.warning(warnstr)
                    # Go to the next line
                    continue

                # Update the FASTA sequence
                tmp_fasta += aa

                # Update the FASTA-like DSSP sequence
                tmp_seqdssp += self.mapping[line[16]]

                # Append a dictionary of data found on this line
                # to the final list
                data.append(\
                    {"#" : int(line[1:5].lstrip(" ")), 
                     "RESIDUE" : int(residue), 
                     "CHAIN" : chain,
                     "AA" : aa,
                     "STRUCTURE" : line[16],
                     "3-TURNS_HELIX" : line[18],
                     "4-TURNS_HELIX" : line[19],
                     "5-TURNS_HELIX" : line[20],
                     "GEOMETRICAL_BEND" : line[21],
                     "CHIRALITY" : line[22],
                     "BETA_BRIDGE_LABEL" : line[23:25],
                     "BP1" : int(line[25:29].lstrip(" ")),
                     "BP2" : line[29:34].lstrip(" ").rstrip(" "),
                     "ACC" : int(line[34:38].lstrip(" ")),
                     "N-H-->O_1" : line[38:50].lstrip(" ").rstrip(" "),
                     "O-->H-N_1" : line[50:61].lstrip(" ").rstrip(" "),
                     "N-H-->O_2" : line[62:73].lstrip(" ").rstrip(" "),
                     "O-->H-N_2" : line[74:85].lstrip(" ").rstrip(" "),
                     "TCO" : float(line[85:91].lstrip(" ")),
                     "KAPPA" : float(line[91:97].lstrip(" ")),
                     "ALPHA" : float(line[97:103].lstrip(" ")),
                     "PHI" : float(line[103:109].lstrip(" ")),
                     "PSI" : float(line[109:115].lstrip(" ")),
                     "X-CA" : float(line[115:122].lstrip(" ")),
                     "Y-CA" : float(line[122:129].lstrip(" ")),
                     "Z-CA" : float(line[129:137].lstrip(" "))})

                # Update the holders of the previous chain
                # and previous residue
                prev_chain = chain
                prev_res = residue

            # Store the FASTA sequence and FASTA-like DSSP
            # sequence for the current chain
            fasta[chain] = tmp_fasta
            seqdssp[chain] = tmp_seqdssp

            # Create the data frame
            df = pd.DataFrame(data).set_index("#").astype(self.DTYPES)

            # Return the data frame and the two dictionaries
            return df, fasta, seqdssp


def run_dssp(dssp_exec,
             pdb_file,
             output_file):
    """Run DSSP.
    """

    # Launch the process
    subprocess.run([dssp_exec,
                    "-i", pdb_file,
                    "-o", output_file])


def get_alphafold_data(uniprot_id,
                       out_pdb_file,
                       out_pae_file):
    """Given a UniProt ID, download data associated with it
    in the AlphaFold database.
    """

    # Generic URL address where AlphaFold PDB models are stored
    ALPHAFOLD_PDB_URL = \
        "https://alphafold.ebi.ac.uk/files/AF-{:s}-F1-model_v4.pdb"

    # Generic URL address where AlphaFold predicted aligned error
    # files are stored
    ALPHAFOLD_PAE_URL = \
        "https://alphafold.ebi.ac.uk/files/" \
        "AF-{:s}-F1-predicted_aligned_error_v4.json"

    # Get the data
    data_pdb = \
        requests.get(ALPHAFOLD_PDB_URL.format(uniprot_id)).text
    data_pae = \
        requests.get(ALPHAFOLD_PAE_URL.format(uniprot_id)).text

    # Write the data to a PDB and a JSON file
    with open(out_pdb_file, "w") as out_pdb, \
         open(out_pae_file, "w") as out_pae:
        out_pdb.write(data_pdb)
        out_pae.write(data_pae)


def get_plddt_per_residue(pdb_file):
    """Get per-residue pLDDT scores from the AlphaFold PDB file.
    """

    # Create an empty list to store the data for each
    # residue
    data = []

    # Get the structure
    structure = bpdb.PDBParser().get_structure("struct", pdb_file)

    # For each chain in the first model (only one model in
    # AlphaFold structures)
    for idx, residue in enumerate(structure[0].get_residues()):

        # Get the residue name
        resname = residue.get_resname()

        # Get only the first atom of the current residue
        # (we are only interest in residue-level properties)
        first_atom = list(residue.get_atoms())[0]

        # Get properties of interest from the ID of the first atom
        struct, mod, chain, (hetflag, resnum, icode), _ = \
            first_atom.get_full_id()
        
        # Get the B-factor for the first atom (= pLDDT score of
        # the residue)
        bfactor = first_atom.get_bfactor()

        # Add a dictionary containing the data for the current
        # residue to the final list
        data.append({"chain" : chain,
                     "resnum" : resnum,
                     "resname" : resname,
                     "pLDDT" : bfactor})

    # Create a data frame with the per-residue data and return it
    return pd.DataFrame(data)


def get_regions_by_plddt(plddt_df, plddt_cutoff):
    """Separate the protein regions having a pLDDT score greater than
    or equal to a pre-defined cut-off from those having a pLDDT score
    lower than the cut-off, and report them in ranges.
    """

    # Define a helper function to define ranges
    def get_ranges(i):
        ranges = []
        for a, b in itertools.groupby(enumerate(i), \
                                      lambda pair: pair[1] - pair[0]):
            b = list(b)
            ranges.append((b[0][1], b[-1][1]))
        return ranges

    # Keep only those residues with pLDDT score >= cut-off
    plddt_greaterequal = \
        plddt_df[plddt_df["pLDDT"] >= plddt_cutoff]["resnum"].tolist()

    # Get the protein regions covered by those residues as ranges
    plddt_greaterequal_ranges = get_ranges(plddt_greaterequal)

    # Convert the ranges whose beginning and end coincide (region
    # with only one residue) to length-1 tuples
    plddt_greaterequal_ranges = \
        [f"{b}-{e}" if b != e else str(b) for (b, e) \
         in plddt_greaterequal_ranges]

    # Keep only those residues with pLDDT score < cut-off
    plddt_lower = \
        plddt_df[plddt_df["pLDDT"] < plddt_cutoff]["resnum"].tolist()

    # Get the protein regions covered by those residues as ranges
    plddt_lower_ranges = get_ranges(plddt_lower)

    # Convert the ranges whose beginning and end coincide (region
    # with only one residue) to length-1 tuples
    plddt_lower_ranges = \
        [f"{b}-{e}" if b != e else str(b) for (b, e) \
         in plddt_lower_ranges]

    # Return the ranges
    return plddt_greaterequal_ranges, plddt_lower_ranges


def main(config,
         wd):

    # Set an empty list to store data for the final data frame
    data = []

    # Set the DSSP output parser
    dssp_parser = DSSPParser()

    # Get the cut-off for the pLDDT score
    plddt_cutoff = config["plddt_cutoff"]
    
    # For each UniProt ID
    for uniprot_id, protein_data in config["uniprot_ids"].items():


        #----------------------- AlphaFold DB ------------------------#


        # Create the directory to store the data
        # about the protein associated with the
        # UniProt ID
        protein_dir_name = protein_data["dir_name"]
        protein_dir = os.path.join(wd, protein_dir_name)
        os.makedirs(protein_dir, exist_ok = True)

        # Download the AlphaFold data into the directory
        out_pdb_file = os.path.join(protein_dir, f"{uniprot_id}.pdb")
        out_pae_file = os.path.join(protein_dir, f"{uniprot_id}.json")
        get_alphafold_data(uniprot_id = uniprot_id,
                           out_pdb_file = out_pdb_file,
                           out_pae_file = out_pae_file)


        #--------------------------- DSSP ----------------------------#


        # Run DSSP on the AlphaFold PDB file
        out_dssp_file = \
            os.path.join(protein_dir, f"{uniprot_id}_dssp.out")
        run_dssp(dssp_exec = config["dssp_exec"],
                 pdb_file = out_pdb_file,
                 output_file = out_dssp_file)

        # Parse the DSSP output
        df, fasta, seqdssp = dssp_parser.parse(out_dssp_file)

        # Select only the columns of interest from the DSSP data frame
        df = df[["CHAIN", "RESIDUE", "AA", "STRUCTURE"]]

        # Change the column names
        df.columns = \
            ["chain", "resnum", "resname_truncated", "secstruc"]

        # Substitute all whitespaces (= coil conformation) with a dash
        # in the column storing the secondary structure conformation
        df = df.replace(to_replace = " ",
                        value = {"secstruc" : "-"})


        #--------------------------- pLDDT ---------------------------#


        # Get the per-residue pLDDT score
        plddt_df = get_plddt_per_residue(out_pdb_file)

        # Get the protein regions with pLDDT score greater than/equal
        # to or lower than the pre-defined pLDDT cut-off
        plddt_greaterequal_ranges, plddt_lower_ranges = \
            get_regions_by_plddt(plddt_df, plddt_cutoff)

        # Append the residue ranges representing these regions
        # to the final list
        data.append(\
                {"uniprot_id" : uniprot_id,
                 "dir_name" : protein_dir_name,
                 f"regions_pLDDT_greaterequal_{plddt_cutoff}" : \
                    ";".join(plddt_greaterequal_ranges),
                 f"regions_pLDDT_lower_{plddt_cutoff}" : \
                    ";".join(plddt_lower_ranges)})


        #------------------ Per-residue data frame -------------------#


        # Create a data frame merging per-residue information about
        # the secodnary structure assigned by DSSP and the pLDDT score
        merged_df = pd.merge(plddt_df, df, on = ["chain", "resnum"])
        merged_df = merged_df.drop("resname_truncated", axis = 1)

        # Save the new data frame to a CSV file
        merged_csv = \
            os.path.join(protein_dir, f"{uniprot_id}.csv")
        merged_df.to_csv(merged_csv, sep = ",", index = False)


    #------------------- Protein regions data frame ------------------#


    # Generate a data frame containing the data about
    # the protein regions with pLDDT scores greater than/equal to
    # or lower than the pre-definied pLDDT cut-off for all proteins
    regions_df = pd.DataFrame(data)

    # Save the data frame in a CSV file
    regions_csv_file = \
        os.path.join(wd, "regions_pLDDT.csv")
    regions_df.to_csv(regions_csv_file, sep = ",", index = False)



if __name__ == "__main__":


    import argparse


    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add the arguments
    c_helpstr = "Configuration file for the curation."
    parser.add_argument("-c", "--config-file",
                        type = str,
                        required = True,
                        help = c_helpstr)

    d_helpstr = "Path to the working directory."
    parser.add_argument("-d", "--work-dir",
                        type = str,
                        default = os.getcwd(),
                        help = d_helpstr)

    # Parse the arguments
    args = parser.parse_args()

    # Get the working directory
    wd = args.work_dir

    # Load the configuration file
    config = yaml.safe_load(open(args.config_file, "r"))

    # Get and save AlphaFold data
    main(config = config,
         wd = wd)

