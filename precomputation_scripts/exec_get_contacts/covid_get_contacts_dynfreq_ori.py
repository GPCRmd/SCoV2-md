#module load Python

import re
import os
import pickle

import glob
import numpy as np
import argparse as ap
import pandas as pd
from random import sample
from json import loads, dump
from sys import stdout
from shutil import copyfile,copyfileobj
import time
#import MDAnalysis as md
#from htmd.ui import *

parser = ap.ArgumentParser(description="this calculates interaction frequencies for given simulation")
parser.add_argument(
   '--dyn',
    dest='dynamics_id',
    nargs='*',
    action='store',
    type=int,
    default=False,
    help='Specify dynamics id(s) for which the matrix will be precomputed. '
)
parser.add_argument(
    '--overwrite',
    dest='repeat_dynamics',
    action='store_true',
    default=False,
    help='(bool) repeat get_dynamic_contacts even if there already exists a results file for this files'
)
#parser.add_argument(
#    '--cores',
#    dest='cores',
#    action='store',
#    default=1,
#    help='number of cores to use in get_dynamic_contacts'
#)
parser.add_argument(
   '-i',
    dest='dyn_dict_path',
    default="/gpcr/users/mariona/input/covid_dyn_dict.data",
    action='store',
    type=str,
    help='Path to the input dictionary with data on the simulations'
)
parser.add_argument(
   '--chains_path',
    dest='fixedchains_path',
    default="/gpcr/users/mariona/covid/exec_get_contacts/dyn_filename_fixedchains.data",
    action='store',
    type=str,
    help='Path to the input dictionary with data on PDB files generted to get data on missing chain indexes.'
)
#parser.add_argument(
#    '--overwrite',
#    action='store_true',
#    dest='overwrite',
#    default=False,
#    help='Overwrites already generated json files.',
#)
parser.add_argument(
    '--ignore_publication',
    action='store_true',
    dest='ignore_publication',
    default=False,
    help='Consider both published and unpublished dynamics.',
)

def get_ligand_selection(lig_list):
    """
    reate a VMD selection line for this ligand. It will be used in --ligand and --sel from 
    get_dynamic_interactions.py
    """
    if not lig_list:
        return("","")
    ligand_sel = ""

    first=True
    for line in lig_list:
        linetab = line.split(":")
        if first: 
            ligand_sel += str("( chain %s and resname %s and resid %s)" % (linetab[0],linetab[1],linetab[2]))
            first = False
        else:
            ligand_sel += str(" or (chain %s and resname %s and resid %s)" % (linetab[0],linetab[1],linetab[2]))
    ligand_apendsel = " or %s" % (ligand_sel)
    ligand_option = "--ligand \"%s\" " % (ligand_sel)
    return(ligand_apendsel, ligand_option)


def get_contacts_dynfreq(dyn_id,mypdbpath,mytrajpath,mytrajid,merge_dynamics,repeat_dynamics,lig_list,cores=4):
    print("Computing dyn %s" % dyn_id)
    dynname = "dyn%s" %dyn_id
    get_contacts_path = "~/bin/getcontacts/"
    files_basepath = "/protwis/sites/files/Precomputed/covid19/get_contacts_files"
    #pharma_path = os.path.join( basepath , "pharmacophores/")
    files_path = os.path.join(files_basepath,"dynfreq_results/dyn%s" % dyn_id)
    if not os.path.isdir(files_path):
        os.makedirs(files_path)

    #Interaction multi-types dictionary
    multi_itypes = {
        'hb' : "hbbb hbsb hbss hbls hblb", # A general category for HB is required
        'wb' : 'wb lwb', # lwb and lwb2 are not needed. I make a posterior division between ligand-residue and residue-residue interactions
        'wb2':'wb2 lwb2',
    }

    #Ligand information extracting
    print("\n")
    (ligand_sel, ligand_option) = get_ligand_selection(lig_list)
    print(ligand_sel)
    print(ligand_option)
    print(lig_list)
    print("\n")
    #Computing dynamic contacts
    dyn_contacts_file=os.path.join(files_path,"%s-%s_dynamic.tsv" % (dynname,mytrajid))
    dyn_contacts_file_merged=os.path.join(files_path,"%s_dynamic.tsv" % dynname)
    if (not os.path.exists(dyn_contacts_file) and not os.path.exists(dyn_contacts_file_merged)) or repeat_dynamics:
        call_py_str=str("python %sget_dynamic_contacts.py         \
        --topology %s  \
        --trajectory %s       \
        --cores %s \
        --sele \"protein%s\"  \
        --itypes all    " % (get_contacts_path, mypdbpath, mytrajpath, cores, ligand_sel) 
        +ligand_option+
        "--output %s" % (dyn_contacts_file)
        )
        print(call_py_str)
        os.system(call_py_str)
    else:
        print("File already exists. Skipping.")
    # Merge dynamic files of this dyn if necessary, and calculate frequencies from this merged dynamic file
    #if merge_dynamics:


def obtain_prot_chains(pdb_filepath,return_model=False):
    try:
        struct = PDB.PDBParser(QUIET=True).get_structure("PDB", pdb_filepath)
        if len(struct) != 1:
            raise ValueError("PDBParser found != 1 model")
    except ValueError:
        print("Warning - PDB parser couñdn't read the PDB file. Manual read is used, which may lead to errors.")
        chains=[]
        fpdb=open(pdb_filepath,'r')
        for line in fpdb:
            if useline2(line):
                chainname=line[21]
                if chainname not in chains:
                    chains.append(chainname)
        model=None
    else:
        model = struct[0]
        model=fixup_resnames(model)
        chains = list(model.child_dict.keys())
    if return_model:
        return chains,model
    else:
        return chains






args = parser.parse_args()
#overwrite=args.overwrite
#exit_on_error=args.exit_on_error
ignore_publication=args.ignore_publication
repeat_dynamics=args.repeat_dynamics
dyn_dict_path=args.dyn_dict_path
dynamics_id=args.dynamics_id
fixedchains_path=args.fixedchains_path


if not os.path.isfile(dyn_dict_path):
    raise ValueError("Input file not found" )
with open(dyn_dict_path, 'rb') as filehandle:  
    dyn_dict = pickle.load(filehandle)
if os.path.isfile(fixedchains_path):
    with open(fixedchains_path, 'rb') as filehandle:  
        fixedchains_d = pickle.load(filehandle)
else:
    print("Warning - file with paths to PDBs with fixed chains not found.")

if dynamics_id:
    dyn_li=[int(d) for d in dynamics_id]
    dyn_dict={k:v for k,v in dyn_dict.items() if k in dyn_li}
if not ignore_publication:
    dyn_dict={k:v for k,v in dyn_dict.items() if v["is_published"]}
if not dyn_dict:
    ValueError("No dynamics found with specified conditions.")
for dyn in sorted(dyn_dict.values(),key=lambda x:x["dyn_id"]):
    dyn_id=dyn["dyn_id"]
    dynfiles=dyn["files"]
    if dyn_id in fixedchains_d:
        mypdbpath=fixedchains_d[dyn_id]
    else:
        mypdbpath=dynfiles["pdb"][0]["path"]
    traj_files=dynfiles["traj"]
    last_traj=len(traj_files)-1
    lig_list=dyn["lig_li_details"]
    merge_dynamics=False
    for idx,mytrajd in enumerate(traj_files):
        if idx==last_traj:
            merge_dynamics=True
        mytrajpath=mytrajd["path"]
        mytrajid=mytrajd["id"]
        get_contacts_dynfreq(dyn_id,mypdbpath,mytrajpath,mytrajid,merge_dynamics,repeat_dynamics,lig_list)
