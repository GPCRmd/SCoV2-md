import re
import os
import pickle
import argparse as ap
from htmd.ui import *


parser = ap.ArgumentParser(description="Checks if some of the input PDB files are missing chain ids and, if its the case, creates new PDB files that have this info.")
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
   '-i',
    dest='dyn_dict_path',
    default="/gpcr/users/mariona/input/covid_dyn_dict.data",
    action='store',
    type=str,
    help='Path to the input dictionary with data on the simulations'
)

parser.add_argument(
    '--ignore_publication',
    action='store_true',
    dest='ignore_publication',
    default=False,
    help='Consider both published and unpublished dynamics.',
)

def modify_filename(mypdbpath):
    filename_ind_start=mypdbpath.rfind("/")
    filename=mypdbpath[filename_ind_start+1:]
    mypath=mypdbpath[:filename_ind_start+1]
    new_filename=os.path.join(mypath,"tmp_chainfix_%s" % filename)
    return new_filename

def fix_fusedchains_problem(mypdbpath):
    mol = Molecule(mypdbpath,validateElements=False)
    chain_li=np.unique(mol.get("chain"))
    chainname_fixed=False
    if len(chain_li)==1 and (chain_li[0] =="" or chain_li[0]==" " ):
        mol = autoSegment(mol,'protein',field="chain", mode="alphabetic")
        chain_li=np.unique(mol.get("chain"))        
        if not (len(chain_li)==1 and (chain_li[0] =="" or chain_li[0]==" " or chain_li[0]=="0" or chain_li[0]==0)):
            chainname_fixed=True
    if chainname_fixed:
        tmpstruc_filepath=modify_filename(mypdbpath)
        mol.write(tmpstruc_filepath)
        return tmpstruc_filepath
    else:
        return False

args = parser.parse_args()
dyn_dict_path=args.dyn_dict_path
dynamics_id=args.dynamics_id
ignore_publication=args.ignore_publication
if not os.path.isfile(dyn_dict_path):
    raise ValueError("Input file not found" )

with open(dyn_dict_path, 'rb') as filehandle:  
    dyn_dict = pickle.load(filehandle)
if dynamics_id:
    dyn_li=[int(d) for d in dynamics_id]
    dyn_dict={k:v for k,v in dyn_dict.items() if k in dyn_li}
if not ignore_publication:
    dyn_dict={k:v for k,v in dyn_dict.items() if v["is_published"]}
if not dyn_dict:
    ValueError("No dynamics found with specified conditions.")

dyn_filename_d={}

for dyn in sorted(dyn_dict.values(),key=lambda x:x["dyn_id"]):
    dyn_id=dyn["dyn_id"]
    dynfiles=dyn["files"]
    mypdbpath=dynfiles["pdb"][0]["path"]
    try:
        mod_filepath=fix_fusedchains_problem(mypdbpath)
    except:
        print("Error in dyn %s" % dyn_id)
        mod_filepath=False
    if mod_filepath:
        dyn_filename_d[dyn_id]=mod_filepath
        print("New file created for dyn %s" % dyn_id)
print("Done!")
out_path="/gpcr/users/mariona/covid/exec_get_contacts/dyn_filename_fixedchains.data"
with open(out_path, "wb") as filehandle:  
    pickle.dump(dyn_filename_d, filehandle)