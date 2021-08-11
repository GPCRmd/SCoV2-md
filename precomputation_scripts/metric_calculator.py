import argparse
import pandas as pd
import numpy as np
from glob import glob
from Bio.Align import substitution_matrices
from Bio import PDB
from Bio.SeqUtils import seq1
import mdtraj as mdt
import logging
import io
import os
import socket
import pickle
import copy
import math
import statistics
import datetime as dt


parser = argparse.ArgumentParser(description="""
Precomputation for variant analyses""") 
parser.add_argument(
    '--debug',
    action='store_true',
    dest='debug',
    default=False,
    help='Debug mode.',
)
parser.add_argument(
    '--overwrite',
    action='store_true',
    dest='overwrite',
    default=False,
    help='Overwrite all results.',
)
parser.add_argument(
   '-i',
    dest='dyn_dict_path',
    default="/gpcr/users/mariona/input/covid_dyn_dict.data",
    action='store',
    type=str,
    help='Path to the input dictionary wit data on the simulations, extracted from the covid platform'
)
parser.add_argument(
   '--matrix_path',
    dest='score_matrices_path',
    default="/home/mariona/Projects/covid19/covid19-genome-analysis/metric/",
    action='store',
    type=str,
    help='Path to the matrices of physico-chemical aminoacid differences.'
)
parser.add_argument(
   '--stride_to_max',
    dest='stride_to_max',
    default=True,
    action='store_false',
    help='If set to true, the simulations are strided until delta=100ps. If this is more than 5000 frames, reduce to delta=10ps.'
)
parser.add_argument(
   '--dyn',
    dest='specify_dyn',
    nargs='*',
    action='store',
    type=int,
    default=False,
    help='Specify dynamics id(s) for which the matrix will be precomputed. '
)
# https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html?highlight=aminoacid
# https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

# model = structure[0]
# chain = model["A"]
# residue = chain[100]
# atom = residue["CA"]


# import ipdb; ipdb.set_trace()

options = parser.parse_args()
debug = options.debug

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
logger = logging.getLogger(__name__)

hostname=socket.gethostname()
if 'MacBook-di-Toni' in hostname:
    out_dir = basedir= dynamics_dir = dynfreq_dir = "scratch"
else:
    if "faramir" in hostname:
        basedir="/home/mariona/gpcrmd_vagrant/shared/sites/files/"
    else:
        basedir="/protwis/sites/files/"
    dynamics_dir = os.path.join(basedir,"Covid19Dynamics")
    dynfreq_dir = os.path.join(basedir,"Precomputed/covid19/get_contacts_files/dynfreq_results")
    out_dir = os.path.join(basedir,"Precomputed/covid19/variant_impact/")
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

def sample_entropy(x, base=2.0, period=2 * np.pi):
    last_ds = (period - np.max(x)) + np.min(x)
    ds = np.diff(np.sort(x))
    ds = np.append(ds, last_ds)
    ds = ds[ds != 0]
    di = 2 * np.pi / ds
    S = -di * np.log(di) / np.log(base)
    return np.sum(S)


class Mutation:
    """Convenience object to represent a mutation string, like "P33I".

    Attributes:
       old: the original AA, 1-letter code
       new: the mutated AA, 1-letter code
       pos: the position
       pdb_aa: can be indicated for the cases where the AA in the PDB that position is not the same as the wt
    """

    def __init__(self, XnnnY,pdb_chain=None,pdb_aa=None,pdb_ca_atom=None):
        self.XnnnY = XnnnY
        wl = len(XnnnY)
        self.old = XnnnY[0]
        self.new = XnnnY[wl-1]
        self.pos = int(XnnnY[1:wl-1])
        self.pdb_ca_atom=pdb_ca_atom
        self.pdb_aa = pdb_aa
        self.pdb_chain = pdb_chain if pdb_chain and pdb_chain.strip() else None

    def __repr__(self):
        return self.XnnnY

def get_aa_charge(myaa):
    positive_aa=["D","E"]
    negative_aa=["K","R","H"]
    charge=0
    if myaa in positive_aa:
        charge=1
    elif myaa in negative_aa:
        charge=-1
    return charge

class MutationEstimator:
    """Main object to compute the impact of point mutations.

    Methods ending in _t return a vector of measurements, one per
    frame (possibly AA-specific, as selected with the "focus"
    method). Otherwise, only one per-trajectory value is returned.

    Methods ending in _m return a mutation-specific measurement,
    i.e. affected by the mutation specified in focus (e.g. PAM
    substitution score). Otherwise, the result is unaffected by
    the mutation (but can still be AA-specific).

    Attributes:
      model, Bio.PDB.Model.Model: the biopython model
      chains, list: the list of chains
      sequences, dict: the sequence of each chain (may be incomplete)
      mdt_traj, mdtraj.core.trajectory.Trajectory: the MDTraj trajectory

    """

    def __init__(self, **kwargs):#
        """Load the trajectory and topology.

        Args:
           traj_id: integer, if set all the other are retrieved from
                    hardcoded paths in covid19 platform, ignoring the following
           pdb: path to the PDB file
           traj: path to the XTC or DCD file
           tsv: path to the .tsv[.xz] precomputed by getContacts
        """

        if "traj_id" in kwargs:
            traj_id = kwargs["traj_id"]
            pdb_l = glob(dynamics_dir+f"/*_dyn_{traj_id}.pdb")
            traj_l = glob(dynamics_dir+f"/*_trj_{traj_id}.*")  # xtc, dcd
            tsv_l = glob(dynfreq_dir+f"/dyn{traj_id}/*.tsv*")  # .xz ?
            if len(pdb_l) != 1 or len(traj_l) != 1 or \
               len(tsv_l) != 1:
                raise FileNotFoundError(
                    f"Inconsistent globbing for trajectory id {traj_id}")
            self.pdb = pdb_l[0]
            self.traj = traj_l[0]
            self.tsv = tsv_l[0]
        else:
            self.pdb = kwargs["pdb"]
            self.traj = kwargs["traj"]
            self.tsv = kwargs["tsv"]
        stride=1
        if "stride" in kwargs:
            stride = kwargs["stride"]
        self.stride=stride
        # Biopython PDB
        struct = PDB.PDBParser(QUIET=True).get_structure("PDB", self.pdb)
        if len(struct) != 1:
            raise ValueError("PDBParser found != 1 model")
        self.model = struct[0]
        self._fixup_resnames()

        # Chains
        self.chains = list(self.model.child_dict.keys())
        if len(self.chains) != 1:
            logger.warning(
                "Multiple chains found, using first, which may not work")

        # Sequences
        ppb = PDB.PPBuilder()
        sequences={}
        for c in self.chains:
            mod_peptides=ppb.build_peptides(self.model[c])
            if mod_peptides:
                sequences[c]=mod_peptides[0].get_sequence()
        self.sequences = sequences
        #self.sequences = {c: ppb.build_peptides(self.model[c])[0].get_sequence()
        #                  for c in self.chains}

        # MDTraj
        self.mdt_traj = mdt.load(self.traj, top=self.pdb, stride=stride)

        # GetContacts
        self._loadContacts()

        # Time-de parameters
        self.sasa_all=None
        self.chi1_all_atomlist=None
        self.chi1_all_angles=None

    def _fixup_resnames(self):
        """Revert non-standard resname nomenclature."""
        fixups = {'CYZ': 'CYS',
                  'CYX': 'CYS',
                  'CYM': 'CYS',
                  'HID': 'HIS',
                  'HIE': 'HIS',
                  'HIP': 'HIS',
                  'HSD': 'HIS',
                  'HSE': 'HIS',
                  'HSP': 'HIS'}
        for r in self.model.get_residues():
            if r.resname in fixups:
                r.resname = fixups[r.resname]

    # def _checkMutN(self, mutation_list):
    #     """Ensure that the mutations match the sequence."""
    #     tot_found = 0
    #     for m in mutation_list:
    #         n_found, n_mismatch = self._findMutation1(m)
    #         tot_found += n_found
    #     if tot_found == 0:
    #         raise ValueError(f"Requested mutation {m} did not match the sequence")

    def _findMutation1(self, mutation, chain=None, icode=None):
        """Iterate over all the residues and return the matching chains, if any

        Args:
          mutation: 
          chain:  (Default value = None)
          icode:  (Default value = None)

        Returns:
          A list of chains where the mutation is found.
        """
        mut_orig_resname= mutobj.pdb_aa or mutobj.old

        n_mismatch = 0
        found_chains = []
        for r in self.model.get_residues():
            _, resid, _ = r.get_id()
            resname3 = r.get_resname()
            resname1 = seq1(resname3)
            reschain = r.get_parent().get_id()

            # Specific chain requested and does not match
            if chain and chain != reschain:
                continue
            if mutation.pos == resid and mut_orig_resname == resname1:
                found_chains.append(reschain)
            elif mutation.pos == resid:
                n_mismatch += 1
                logger.warning(
                    f"Mismatched mutation string: {mutation} but resid {mutation.pos} of chain `{reschain}' is {resname3}")
        return found_chains




    def _loadContacts(self):
        """Load the matrix generated by getContacts."""
        nrows = 20000 if debug else None

        chunksize = 10 ** 6

        is_strided=False
        keep_frames=None
        if self.stride and self.stride >1:
            is_strided=True
            total_traj_frames=(self.mdt_traj.n_frames-1)*self.stride
            keep_frames = np.arange(0, int(total_traj_frames), int(stride))

        readdf = pd.read_table(self.tsv,
                           skiprows=2,
                           nrows=nrows,
                           header=None,
                           chunksize=chunksize,
                           names=["frame", "interaction_type",
                                  "atom_1", "atom_2"])
        isFirst=True
        with readdf as reader:
            for chunk in reader:
                if is_strided:
                    chunk=chunk[chunk["frame"].isin(keep_frames)]
                if isFirst:
                    df=chunk
                    isFirst=False
                else:
                    df=df.append(chunk)
        df=df.reset_index()
#        df = pd.read_table(self.tsv,
#                           skiprows=2,
#                           nrows=nrows,
#                           header=None,
#                           names=["frame", "interaction_type",
#                                  "atom_1", "atom_2"])
#        #if stride, remove from the df the info of unneeded frames
#        if self.stride:
#            if self.stride >1:
#                total_traj_frames=(self.mdt_traj.n_frames-1)*self.stride
#                keep_frames = np.arange(0, int(total_traj_frames), int(stride))
#                df=df[df["frame"].isin(keep_frames)]

        # These two are ok but kill machines for memory occupation
        # df[['chain_1','resname_1','resid_1','name_1']] = df.atom_1.str.split(":", expand=True)
        # df[['chain_2','resname_2','resid_2','name_2']] = df.atom_2.str.split(":", expand=True)
        a1 = io.StringIO("\n".join(df.atom_1))
        del df['atom_1']
        df1 = pd.read_table(a1, sep=":", header=None, names=[
                            'chain_1', 'resname_1', 'resid_1', 'name_1'])
        a2 = io.StringIO("\n".join(df.atom_2))
        del df['atom_2']
        df2 = pd.read_table(a2, sep=":", header=None, names=[
                            'chain_2', 'resname_2', 'resid_2', 'name_2'])
        # import ipdb; ipdb.set_trace()
        self.contacts = pd.concat([df, df1, df2], axis=1)

    def _interpolate(self, x):
        """Interpolate a given vector to the length of the trajectory."""
        NS = len(x)
        NF = self.mdt_traj.n_frames
        out = np.interp(np.linspace(0, NS-1, NF), range(NS), x)
        return out

    def focus(self, mut, chain=None, icode=None):
        """Start analyzing a given mutation.

        Args:
          mut: the mutation, as a string or Mutation object.
          chain:  (Default value = None)
          icode:  (Default value = None)

        Returns:

        """
        if type(mut) == str:
            mut = Mutation(mut)

        m_chain_list = self._findMutation1(mut, chain, icode)
        multiple_chains=False
        if len(m_chain_list) > 1:
            if mut.pdb_chain:
                m_chain_list=[e for e in m_chain_list if e.strip()==mut.pdb_chain]
                if len(m_chain_list) > 1:
                    multiple_chains=True
            else:
                multiple_chains=True
        if multiple_chains:
            logger.warning("Several chains matching, not implemented yet")
        m_chain = m_chain_list[0]
        #logger.info("Matched chain "+m_chain)

        try:
            self.m_biopython_residue = self.model[m_chain][mut.pos]
        except KeyError:
            raise ValueError(f"Residue {mut.pos} not in PDB")

        if mut.pdb_ca_atom:
            self.m_mdt_residue=self.mdt_traj.topology.atom(mut.pdb_ca_atom).residue
            self.m_mdt_atoms = np.array([a.index for a in self.m_mdt_residue.atoms])
        else:
            self.m_mdt_atoms = self.mdt_traj.topology.select(f"residue {mut.pos}")
            self.m_mdt_residue = self.mdt_traj.topology.atom(
                self.m_mdt_atoms[0]).residue

        found_residues={ self.mdt_traj.topology.atom(atom).residue for atom in  self.m_mdt_atoms}
        if len(found_residues)>1:
            logger.warning("Waning - Multiple residues selected for mutation %s" % mut)
        logger.info(found_residues)

        self.mut = mut
        self.m_chain_list = m_chain_list

    def get_substitution_matrix_score_m(self, matrix="BLOSUM90"):
        """Returns the fixed (AA-to-AA) substitution score.

        The full list of scores, with references, is in the data/README.md file. 

        See also https://biopython.org/docs/dev/api/Bio.Align.substitution_matrices.html,
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html, and
        http://stevegallik.org/cellbiologyolm_Ex02_P03.html . 

        Args:
          matrix:  (Default value = "BLOSUM90")

        Returns:
          The substitution penalty or change in the selected property.
        """
        try:
            # Search in the standard place
            sm = substitution_matrices.load(matrix)
        except FileNotFoundError:
            # Search our local supplied matrices
            #path = os.path.realpath(__file__)
            #directory = os.path.dirname(path)
            directory=score_matrices_path
            datafile = os.path.join(directory, "data", matrix)
            sm = substitution_matrices.read(datafile)
        return sm[self.mut.old, self.mut.new]

    def get_charge_diff(self):
        """
        Returns, in absolute value, the charge difference between the mutated aminoacid and the original
        """
        final_charge=get_aa_charge(me.mut.new) - get_aa_charge(me.mut.old)
        return final_charge

    def get_hydro_diff(self):
        """
        Returns, in absolute value, the hydrophobicity difference between the mutated aminoacid and the original, based on the Kyte Doolittle scale.
        """
        hydrophobicity_Kyte_Doolittle={'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2};
        return hydrophobicity_Kyte_Doolittle[me.mut.new] - hydrophobicity_Kyte_Doolittle[me.mut.old]
        

    def _compute_sasa_allres(self, stride=None):
        """Compute the Shrake and Rupley solvent-accessible surface. It then assigns to self the SASA (nm^2) as a function of time.   

        See https://mdtraj.org/1.9.4/examples/solvent-accessible-surface-area.html

        Shrake A, Rupley JA. Environment and exposure to solvent of
        protein atoms. Lysozyme and insulin. Journal of Molecular
        Biology. 1973 Sep 15;79(2):351–71. 10.1016/0022-2836(73)90011-9

        Args:
          stride, int: frame subsetting.
   
        """
        logger.info("Computing SASA!! :D ")
        if not stride:
            stride = 1000 if debug else 1

        sasa_all = mdt.shrake_rupley(self.mdt_traj[::stride], mode="residue")
        self.sasa_all=sasa_all



    def get_sasa_t(self, stride=None):
        """Compute the Shrake and Rupley solvent-accessible surface.

        See https://mdtraj.org/1.9.4/examples/solvent-accessible-surface-area.html

        Shrake A, Rupley JA. Environment and exposure to solvent of
        protein atoms. Lysozyme and insulin. Journal of Molecular
        Biology. 1973 Sep 15;79(2):351–71. 10.1016/0022-2836(73)90011-9

        Args:
          stride, int: frame subsetting.

        Returns:
          The SASA (nm^2)  of the mutated (focus) residue as a function of time .If stride > 1, missing frames 
          are interpolated.
        """
        if self.sasa_all is None:
            self._compute_sasa_allres()

        sasa_all=self.sasa_all
        sasa_res = sasa_all[:, self.m_mdt_residue.index]
        return self._interpolate(sasa_res)

    def get_contacts_t(self, mode="vdw"):
        """Return the precomputed contacts of given type.

        Args:
          mode, str: one of those supported: vdw, sb, pc, ps, ts, hbbb, hbsb, hbss

        Returns:
          the number of contacts as a function of time. For each frame,
          return the number of contacts marked. Note that bigger residues
          make more contacts.
        """
        sel = self.mut.pos
        if (mode=="all"):
            dat = self.contacts
        elif (mode=="hb"):
            dat = self.contacts.loc[(self.contacts.interaction_type=="hbbb") | (self.contacts.interaction_type=="hbsb") | (self.contacts.interaction_type=="hbss") | (self.contacts.interaction_type=="hb")]
        else:
            dat = self.contacts.loc[self.contacts.interaction_type == mode]
        rdat = dat[(dat.resid_1 == sel) | (dat.resid_2 == sel)]
#        if (mode=="hb" or mode =="all"):
        rdat = rdat.drop_duplicates(subset=['frame', 'chain_1', 'resid_1', 'chain_2', 'resid_2'], keep='first')
        frc = rdat.groupby("frame").count()['interaction_type']
        if self.stride and self.stride>1:
            frame_range=range(0,int(self.mdt_traj.n_frames)*int(stride),int(self.stride))
        else:
            frame_range=range(self.mdt_traj.n_frames)
        frc0 = frc.reindex(frame_range, fill_value=0)
        return frc0.to_numpy()

    def get_rmsd_aa_t(self):
        """Return the RMSD of the aminoacid w.r.t. the first frame"""

        # TODO: align?
        rmsd = mdt.rmsd(self.mdt_traj,
                        self.mdt_traj,
                        atom_indices=self.m_mdt_atoms)
        return rmsd

    def get_rmsf_aa_t(self):
        """Return the RMSD of the aminoacid w.r.t. the first frame"""

#        rmsf = mdt.rmsf(self.mdt_traj,
#                        self.mdt_traj,
#                        atom_indices=self.m_mdt_atoms)

        #We 1st align selected atoms
        superposed_traj=copy.deepcopy(self.mdt_traj)
        superposed_traj.superpose(atom_indices=self.m_mdt_atoms,reference=superposed_traj)
        rmsf = mdt.rmsf(target=superposed_traj,
                        reference=None,
                        atom_indices=me.m_mdt_atoms)
        del superposed_traj


        return rmsf

    def _compute_chi1_allres(self):
        """Computes the chi1 angle of all aminoacids."""
        alists, chi1s = mdt.compute_chi1(self.mdt_traj)
        self.chi1_all_atomlist=alists
        self.chi1_all_angles=chi1s


    def get_chi1_t(self):
        """Return the chi1 angle of the aminoacid.

        Returns False if there is no such chi angle."""

        if self.chi1_all_atomlist is None or self.chi1_all_angles is None:
            self._compute_chi1_allres()
        alists=self.chi1_all_atomlist
        chi1s=self.chi1_all_angles

        focus_atoms = set(self.m_mdt_atoms)
        # Pick the correct AA
        for i, al in enumerate(alists):
            if set(al).issubset(focus_atoms):
                return chi1s[:, i]
        logger.warning("Focused atom does not have chi1 angle")
        return False

    def get_rotameric_entropy(self, base=2.0, period=2*np.pi):
        """Return the 'rotameric chi1 entropy' of the aminoacid.

        This is my definition (Toni's) of a sample entropy and may be wrong.

        Args:
          base: float, the base for the entropy calculation (def. bits)
          period: the periodicity (default 2 pi)
        """
        chi1 = self.get_chi1_t()
        return sample_entropy(chi1, base, period)

def get_or_create_analysis_res(out_dir,analysis,dyn_id,traj_id,overwrite,result_is_pandas=False,include_traj=False):
    analysis_dir=os.path.join(out_dir,analysis)
    if not os.path.isdir(analysis_dir):
        os.mkdir(analysis_dir)
    if result_is_pandas:
        ext="csv"
        my_result=pd.DataFrame()
    else:
        ext="data"
        my_result={}
    if include_traj:
        out_file=os.path.join(analysis_dir,"dyn_%s_traj_%s.%s" % (dyn_id,traj_id,ext))
    else:
        out_file=os.path.join(analysis_dir,"dyn_%s.%s" % (dyn_id,ext))
    if not overwrite and os.path.isfile(out_file):
        if result_is_pandas:
            my_result=pd.read_csv(out_file)
        else:
            with open(out_file, 'rb') as filehandle:     
                my_result=pickle.load(filehandle)
    return(my_result,out_file)

def transform_res(res_array,stride,head=False):
    if not stride:
        stride=1
    array_p=[[i*stride,float(val)] for i,val in enumerate(res_array)]
    if head:
        array_p.insert(0,head) 
    return array_p

def calculate_stride(delta,framenum):
    stride=delta*10 #the desired delta is 100ps (0.1 ns)
    final_num_frames=math.floor(framenum/stride)
    if final_num_frames> 5000:
        stride=delta*100 #the desired delta is 10ps (0.01 ns)
    stride=math.ceil(stride)
    return stride

def compute_variables(me,analysis,stride):

    if analysis.startswith("contact_"):
        contact_type=analysis.split("_")[1]
        myresult=me.get_contacts_t(contact_type)
        myresult_p=transform_res(myresult,stride)
    elif analysis == "rmsd":
        myresult=me.get_rmsd_aa_t()
        myresult_p=transform_res(myresult,stride)
    elif analysis == "sasa":
        myresult=me.get_sasa_t()
        myresult_p=transform_res(myresult,stride)
    elif analysis == "chi1":
        myresult=me.get_chi1_t()
        myresult_p=transform_res(myresult,stride)
    elif analysis == "rmsf":
        myresult=me.get_rmsf_aa_t()
        atom_li=[" ".join([str(a.index),a.name]) for a in me.m_mdt_residue.atoms]
        myresult_p=[[atom_li[i],float(val)] for i,val in enumerate(myresult)]
    return myresult_p

def average(lst): 
    return sum(lst) / len(lst) 

def get_time():
    mytime=dt.datetime.today()
    myhour=mytime.hour
    mymin=mytime.minute
    return (myhour,mymin)


if __name__ == "__main__":

    #me = MutationEstimator(traj_id=5)
    #P33I = Mutation("P33I")
    #me.focus(P33I)

    overwrite=options.overwrite
    dyn_dict_path=options.dyn_dict_path
    score_matrices_path=options.score_matrices_path
    stride_to_max=options.stride_to_max
    specify_dyn=options.specify_dyn


    all_analysis_frame=["rmsd","sasa","chi1",'contact_hb', 'contact_sb', 'contact_hp', 'contact_pc', 'contact_ps', 'contact_ts', 'contact_vdw', 'contact_wb', 'contact_wb2', 'contact_all']
    all_analysis_other=["rmsf"]

    if not os.path.isfile(dyn_dict_path):
        raise Exception("Database data not found")
    with open(dyn_dict_path, 'rb') as filehandle:  
        dyn_dict = pickle.load(filehandle)

    if specify_dyn:
        dyn_li=[int(d) for d in specify_dyn]
        dyn_dict={k:v for k,v in dyn_dict.items() if k in dyn_li}

    for dyn_id,dyn in sorted(dyn_dict.items(),key=lambda x:x[0]):
        sim_files=dyn["files"]
        pdb_files=sim_files["pdb"]
        delta=dyn["delta"]
        if not delta:
            if int(dyn_id)==5:
                delta=0.2
            else:
                logger.info("No delta")
                continue
        if not dyn["is_published"]:
            continue
        if not pdb_files:
            continue
        pdb_path=pdb_files[0]["path"]
        traj_files=sim_files["traj"]
        traj_done=set()
        for trajfile in traj_files:
            (hour_start_traj,min_start_traj)=get_time()

            traj_id=trajfile["id"]
            if traj_id in traj_done:
                continue
            traj_done.add(traj_id)
            logger.info("\nDyn %s - Traj %s" % (dyn_id,traj_id))
            framenum=trajfile["framenum"]
            stride=False
            if stride_to_max:
                stride=calculate_stride(delta,framenum)
            traj_path=trajfile["path"]
            contacts_path=os.path.join(dynfreq_dir,"dyn%(dyn_id)s/dyn%(dyn_id)s-%(traj_id)s_dynamic.tsv"%{"dyn_id":dyn_id,"traj_id":traj_id})
            try:
                me = MutationEstimator(pdb=pdb_path,traj=traj_path,tsv=contacts_path,stride=stride)
            except Exception as e:
                logger.info("Skipping %s - %s"% (dyn_id,e))
                continue
            dyn_results_impact={}
            (mydyn_res_mut_impact,mut_impact_out_file)=get_or_create_analysis_res(out_dir,"mut_impact",dyn_id,traj_id,overwrite)
            (mydyn_var_summary,summary_out_file)=get_or_create_analysis_res(out_dir,"summary",dyn_id,traj_id,overwrite,include_traj=True)
            dyn_results_frame={}
            for analysis in all_analysis_frame:
                (my_result,out_file)=get_or_create_analysis_res(out_dir,analysis,dyn_id,traj_id,overwrite,True,True)
                dyn_results_frame[analysis]={}
                dyn_results_frame[analysis]["filepath"]=out_file
                dyn_results_frame[analysis]["result"]=my_result

            for (model_pos,model_chain),pos_data in sorted(dyn['model_to_seq'].items(), key=lambda x:x[0][0]):
                pos_computed=False
                logger.info("Model pos: %s" % model_pos)
                pos_model_ca_atom=pos_data["model_details"]["ca_atom_index"]
                pos_model_chain=pos_data["model_details"]["chain"]
                pos_model_aa=pos_data['model_details']["aa"]
                for protname, prot_pos_data in pos_data["finprot_seq_details"].items():
                    prot_var_data=prot_pos_data["var_data"]
                    if protname not in mydyn_var_summary:
                        mydyn_var_summary[protname]={}
                    if protname not in mydyn_res_mut_impact:
                        mydyn_res_mut_impact[protname]={}
                    for varname, var_data in prot_var_data.items():
                        #logger.info("Variant %s"%varname)
                        fromaa=var_data["resletter_from"]
                        toaa=var_data["resletter_to"]
                        fp_seq_pos=var_data["resid"]
                        mutation_nm_model="%s%s%s"%(fromaa,model_pos,toaa)
                        mutobj=Mutation(mutation_nm_model,pdb_chain=pos_model_chain,pdb_aa=pos_model_aa,pdb_ca_atom=pos_model_ca_atom) 
                        try:
                            me.focus(mutobj)
                        except Exception as e:
                            logger.info("ERROR in variant %s: %s" % (varname,e))
                            continue

                        if fp_seq_pos not in mydyn_res_mut_impact[protname]:
                            mydyn_res_mut_impact[protname][fp_seq_pos]={}
                        mydyn_res_mut_impact[protname][fp_seq_pos][varname]={}
                        blosum90_score=me.get_substitution_matrix_score_m() #maybe I can remove this?
                        mydyn_res_mut_impact[protname][fp_seq_pos][varname]["blosum90"]=blosum90_score #maybe I can remove this?

                        if (fp_seq_pos not in mydyn_var_summary[protname]):
                            mydyn_var_summary[protname][fp_seq_pos]={"time_dep":{}, "variants":{}}

                        mydyn_var_summary[protname][fp_seq_pos]["variants"][varname]={}
                        mut_effect_matrix_params=['BLOSUM90', 'CHARGE', 'EPSTEIN', 'EXEX', 'GRANTHAM', 'MIYATA', 'SNEATH', 'HP_EISENBERG_WEISS', 'HP_ENGLEMAN', 'HP_HESSA', 'HP_HOOP_WOODS', 'HP_JANIN', 'HP_KYTE_DOOLITTLE', 'HP_KYTE_DOOLITTLE_2', 'HP_MOON_FLEMING', 'HP_WIMLEY_WHITE', 'HP_ZHAO_LONDON']
                        for me_mp in mut_effect_matrix_params:
                            me_result=me.get_substitution_matrix_score_m(me_mp)
                            param_ref="variant_"+me_mp.lower()
                            if me_mp.startswith("HP_"):
                                if "hydrophobicity" not in mydyn_var_summary[protname][fp_seq_pos]["variants"][varname]:
                                    mydyn_var_summary[protname][fp_seq_pos]["variants"][varname]["hydrophobicity"]={}
                                mydyn_var_summary[protname][fp_seq_pos]["variants"][varname]["hydrophobicity"][param_ref]= me_result
                            else:
                                mydyn_var_summary[protname][fp_seq_pos]["variants"][varname][param_ref]= me_result


                        if not pos_computed: #computaions considered here are not affected by resulting aa of the variant
                            for analysis,analysis_data in dyn_results_frame.items():
                                my_results=analysis_data["result"]
                                if str(fp_seq_pos) not in my_results.columns:
                                    (hour_start,min_start)=get_time()
                                    try:
                                        analysis_res=compute_variables(me,analysis,stride)
                                    except Exception as e:
                                        logger.warning("Error when computing %s for dyn %s: %s" % (analysis,dyn_id, e))
                                        continue
                                    my_results[str(fp_seq_pos)]=[e[1] for e in analysis_res]
                                    out_file=analysis_data["filepath"]
                                    my_results.to_csv(out_file)
    #                                with open(out_file, 'wb') as filehandle:     
    #                                    pickle.dump(my_results, filehandle)
                                    (hour_end,min_end)=get_time()
                                    total_min=(hour_end - hour_start)*60 + (min_end - min_start)

                                    logger.info("Saved %s (Time: %s min)" % (analysis,total_min))
                                if analysis not in mydyn_var_summary[protname][fp_seq_pos]["time_dep"]:
                                    position_results=my_results[str(fp_seq_pos)].values.tolist()
                                    mydyn_var_summary[protname][fp_seq_pos]["time_dep"][analysis]={}
                                    mydyn_var_summary[protname][fp_seq_pos]["time_dep"][analysis]["average"]=average(position_results)
                                    mydyn_var_summary[protname][fp_seq_pos]["time_dep"][analysis]["sd"]=statistics.pstdev(position_results)

                            for analysis in all_analysis_other:
                                analysis_dir=os.path.join(out_dir,analysis)
                                if not os.path.isdir(analysis_dir):
                                    os.mkdir(analysis_dir)
                                result_filepath=os.path.join(analysis_dir,"dyn_%s_traj_%s_pos_%s_.data" % (dyn_id,traj_id,fp_seq_pos))
                                analysis_res=False
                                if overwrite or not os.path.isfile(result_filepath):
                                    try:
                                        analysis_res=compute_variables(me,analysis,stride)
                                    except Exception as e:
                                        logger.warning("Error when computing %s for dyn %s: %s" % (analysis,dyn_id, e))
                                        continue
                                    with open(result_filepath, 'wb') as filehandle:     
                                        pickle.dump(analysis_res, filehandle)
                                    logger.info("Saved %s" % analysis)
                                if analysis not in mydyn_var_summary[protname][fp_seq_pos]["time_dep"]:
                                    if not analysis_res:
                                        with open(result_filepath, 'rb') as filehandle:     
                                            analysis_res=pickle.load(filehandle)
                                    analysis_values=[e[1] for e in analysis_res]
                                    mydyn_var_summary[protname][fp_seq_pos]["time_dep"][analysis]={}
                                    mydyn_var_summary[protname][fp_seq_pos]["time_dep"][analysis]["average"]=average(analysis_values)
                                    mydyn_var_summary[protname][fp_seq_pos]["time_dep"][analysis]["sd"]=statistics.pstdev(analysis_values)
                                    logger.info("Obtained summary of %s" % analysis)

                            pos_computed=True

                #save mut impact
                with open(mut_impact_out_file, 'wb') as filehandle:     
                    pickle.dump(mydyn_res_mut_impact, filehandle)
                logger.info("Saved mutation impact")

                #save per-variant summary
                with open(summary_out_file, 'wb') as filehandle:     
                    pickle.dump(mydyn_var_summary, filehandle)
                logger.info("Saved per-variat summary")
            (hour_end_traj,min_end_traj)=get_time()
            total_min_traj=(hour_end_traj - hour_start_traj)*60 + (min_end_traj - min_start_traj)
            logger.info("Total time for traj.: %s min" % total_min_traj)
