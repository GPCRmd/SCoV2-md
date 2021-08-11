
import urllib
import copy
import pickle
import os
import mdtraj as md
from Bio import pairwise2
import requests
import gc
from Bio import PDB
from htmd.ui import *

def useline2(line):
    '''returns True if line starts with ATOM, or HETATM with a resname included in the d dictionary''' 
    if line.startswith('ATOM') or line.startswith('HETATM'):
        trykey=line[17:21]
        trykey=trykey.strip()
        if trykey in d.keys():
            return True
        else:
            return False #this heteroatom is not useful

    else:
        return False

def fixup_resnames(mymodel):
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
    for r in mymodel.get_residues():
        if r.resname in fixups:
            r.resname = fixups[r.resname]
    return mymodel

#def obtain_prot_chains_pre(pdb_filepath):
#    chain_name_li=[]
#    fpdb=open(pdb_filepath,'r')
#    for line in fpdb:
#        if useline2(line):
#            chainname=line[21]
#            if chainname not in chain_name_li:
#                chain_name_li.append(chainname)
#    return list(chain_name_li)

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

def obtain_chain_breaks_htmd(mypdb_filepath):
    mol = Molecule(mypdb_filepath)
    chain_li=np.unique(mol.get("chain"))
    chainname_missing=False
    if len(chain_li)==1 and (chain_li[0] =="" or chain_li[0]==" "):
        mol = autoSegment(mol,'protein',field="chain", mode="alphabetic")
        chain_li=np.unique(mol.get("chain"))
        if len(chain_li)==1 and (chain_li[0] =="" or chain_li[0]==" " or chain_li[1]):
            chainname_missing=True


    chain_seq=[]
    chain_count=0
    chain_breaks=[]
    for chainname in chain_li:
        prot_sel="chain %s and name CA" % chainname
        if chainname_missing:
            chainname=chr(ord('A')+chain_count)
            prot_sel="name CA"
        if not chainname:
            continue
        resSeq_li=mol.get("resid" ,sel=prot_sel)
        if len(resSeq_li)==0:
            continue
        chain_breaks.append(resSeq_li[-1]-1)
        chain_count+=1
    chain_breaks=sorted(list(set(chain_breaks)))
    return chain_breaks



def obtain_chain_breaks(struc):
    """
        Returns list of residue index that separate chains. The residue index indicated is the last residue of the chain.
        For that, checks which consecutive residues are far apart (distance between CAs higher than 5A).
    """
    if type(struc)==str:
        struc=md.load(struc)
    all_resindex=list(set([e.index for e in struc.topology.residues if e.name in d]))
    consec_resindex_pairs=[[residx , all_resindex[i+1]] for (i,residx) in enumerate(all_resindex) if i+1 < len(all_resindex)]
    (dists,res_p)=md.compute_contacts(traj=struc,contacts=consec_resindex_pairs,scheme="ca")
    dists_str=dists[0] #we don't have multiple frames
    chain_break_idxs=list()
    for i, dist in enumerate(dists_str):
        if dist > 0.5:
            (residx0,residx1)=res_p[i]
            # We need to validate that the residues found are really consecutive residues that are separated. compute_contacts ignores residues that are not considered to be protein. Thus, some protein residues that have unclear names in the PDB (ex CYZ) are ignored. This may lead to believe that residues separated by this type of prot residues (for ex. CYZ) are consecutive residues that are far apart, which is not the case.
            if residx1 != residx0+1: # some residue between them
                #check if ahe residues between are protein
                between_res=range(residx0+1,residx1)
                allprot_between=True
                for ridx in between_res:
                    if struc.topology.residue(ridx).name not in d: 
                        allprot_between=False
                if allprot_between:
                    continue
            #print(dist,res_p[i],struc.topology.residue(res_p[i][0]),struc.topology.residue(res_p[i][1]))
            chain_break_idxs.append(residx0)
    return chain_break_idxs



def mdtraj_resindex_to_resseq(struc,resindex):
    return struc.topology.atom(struc.topology.select("resid %s" % resindex)[0]).residue.resSeq



def extract_PDB_seq_mdtraj(mypdb_filepath,include_ca_index=False,has_null_chains=False,chain_li=None):
    struc=md.load(mypdb_filepath)
    if has_null_chains:
        chain_break_idxs=obtain_chain_breaks(struc)
        missing_chain_num=0
    protchain_idx=0
    chain_seq=[]
    chainname=False
    for chain in struc.topology.chains:
        tablepdb=[]
        pdb_sequence=""
        for residue in chain.residues:
            resname=residue.name
            if resname in d:  #is protein
                residx=residue.index
                resSeq=residue.resSeq
                res1let=d[resname]
                pdb_sequence+=res1let
                try:
                    ca_index=residue.atom("CA").index
                except KeyError:
                    ca_index=None
                if include_ca_index:
                    tablepdb.append([res1let,int(resSeq),ca_index])
                else:
                    tablepdb.append([res1let,int(resSeq)])
                if has_null_chains:
                    if residx in chain_break_idxs:
                        if pdb_sequence:
                            chainname=chr(ord('A')+missing_chain_num)
                            chain_info={"pdbseq":pdb_sequence,"tablepdb":tablepdb}
                            chain_info["chainname"]=chainname
                            chain_seq.append(chain_info)
                        tablepdb=[]
                        pdb_sequence=""
                        missing_chain_num+=1
        if pdb_sequence:
            chain_info={"pdbseq":pdb_sequence,"tablepdb":tablepdb}
            chainname=None
            if has_null_chains:
                chainname=chr(ord('A')+missing_chain_num)
            else:
                if chain_li:
                    chainname=chain_li[protchain_idx]
                    protchain_idx+=1
            if chainname:
                chain_info["chainname"]=chainname
            chain_seq.append(chain_info)
    return chain_seq

def extract_PDB_seq_biopdb(mypdb_filepath,model,include_ca_index=False,has_null_chains=False):
    if has_null_chains:
        chain_break_idxs=obtain_chain_breaks(mypdb_filepath)
        chain_break_resseq=[mdtraj_resindex_to_resseq(struc,e) for e in chain_break_idxs]
        missing_chain_num=0

    chain_seq=[]
    chainname=False
    for chainname,chain in model.child_dict.items():
        tablepdb=[]
        pdb_sequence=""
        for residue in chain:
            resname=residue.get_resname()
            if resname in d:
                resSeq=residue.get_id()[1]
                if resSeq==0:
                    try:
                        residue=chain[1]
                        resname=residue.get_resname()
                        resSeq=1
                    except KeyError:
                        pass
                res1let=d[resname]
                pdb_sequence+=res1let                    
                ca_index=None
                for atom in residue:
                    if atom.get_name()=="CA":
                        ca_index=atom.get_serial_number() -1
               # print(residue.get_resname(),list(residue.get_id())[1],ca) 
                if include_ca_index:
                    tablepdb.append([res1let,int(resSeq),ca_index])
                else:
                    tablepdb.append([res1let,int(resSeq)])
                if has_null_chains:
                    if resSeq in chain_break_resseq:
                        if pdb_sequence:
                            chainname=chr(ord('A')+missing_chain_num)
                            chain_info={"pdbseq":pdb_sequence,"tablepdb":tablepdb}
                            chain_info["chainname"]=chainname
                            chain_seq.append(chain_info)
                        tablepdb=[]
                        pdb_sequence=""
                        missing_chain_num+=1

        if pdb_sequence:
            if has_null_chains:
                chainname=chr(ord('A')+missing_chain_num)
            chain_info={"pdbseq":pdb_sequence,"tablepdb":tablepdb}
            chain_info["chainname"]=chainname
            chain_seq.append(chain_info)
    return chain_seq




def extract_PDB_seq_htmd(mypdb_filepath,include_ca_index=False):
    mol = Molecule(mypdb_filepath)
    chain_li=np.unique(mol.get("chain"))
    chainname_missing=False
    if len(chain_li)==1 and (chain_li[0] =="" or chain_li[0]==" " or chain_li[1]):
        mol = autoSegment(mol,'protein',field="chain", mode="alphabetic")
        chain_li=np.unique(mol.get("chain"))
        if len(chain_li)==1 and (chain_li[0] =="" or chain_li[0]==" " or chain_li[1]):
            chainname_missing=True


    chain_seq=[]
    chain_count=0
    for chainname in chain_li:
        prot_sel="chain %s and name CA" % chainname
        if chainname_missing:
            chainname=chr(ord('A')+chain_count)
            prot_sel="name CA"
        if not chainname:
            continue
        tablepdb=[]
        pdb_sequence=""

        resSeq_li=mol.get("resid" ,sel=prot_sel)
        resname_li=mol.get("resname" ,sel=prot_sel)
        ca_index_li=mol.get("index" ,sel=prot_sel)
        if len(resSeq_li)==0:
            continue
        for ri,resSeq in enumerate(resSeq_li):
            resname=resname_li[ri]
            if resname in d:
                res1let=d[resname]
                ca_index=ca_index_li[ri]
                pdb_sequence+=res1let                    
                if include_ca_index:
                    tablepdb.append([res1let,int(resSeq),ca_index])
                else:
                    tablepdb.append([res1let,int(resSeq)])
        if pdb_sequence:
            chain_info={"pdbseq":pdb_sequence,"tablepdb":tablepdb,"chainname":chainname}
            chain_seq.append(chain_info)
            chain_count+=1
    return chain_seq




def extract_PDB_seq(mypdb_filepath,include_ca_index=False):
    (chain_li,model)=obtain_prot_chains(mypdb_filepath,True)
    if chain_li == [" "] or chain_li == [""]: #For PDB files without chain name info
        chain_seq = extract_PDB_seq_htmd(mypdb_filepath,include_ca_index=include_ca_index)
        #chain_seq = extract_PDB_seq_mdtraj(mypdb_filepath,include_ca_index=include_ca_index,has_null_chains=True) #
        #chain_seq = extract_PDB_seq_biopdb(mypdb_filepath,model,include_ca_index,True)  #Ibiopython do not detect the two chains (ignores repeated resids)
    else:
        if model:
            chain_seq = extract_PDB_seq_biopdb(mypdb_filepath,model,include_ca_index)
        else:
            chain_seq = extract_PDB_seq_mdtraj(mypdb_filepath,include_ca_index=include_ca_index,chain_li=chain_li) 
    return chain_seq

def checkpdb_ngl(name_of_file,segid,start,stop,chain):
    '''Get sequence from a PDB file in a given interval defined by a combination of Segment Identifier (segid), starting residue number (start), end residue number (stop), chain identifier (chain). All can be left in blank. Returns 1) a list of minilist: each minilist has the resid and the aminoacid code. 2) a string with the sequence.'''
    fpdb=open(name_of_file,'r')
    cpos=0 #current residue position
    ppos=0 #previous residue position
    ppos2='0' #previous position after converting hexadecimals to decimals
    pchain='' #previous chain
    seqplain=list() #list of minilist. each minilist contains the residue number and its aminoacid.
    flag=0
    hexflag=0
    pfields=['','' ,'','AAA','Z','0','0','0','0','']
    for line in fpdb:
        if useline2(line):
            fields=[ '','' ,'' ,line[17:21],line[21],line[22:27],line[31:39],line[39:47],line[47:55],line[72:77]] 
            #fields[3]:Aminoacid code, fields[4]:chain, fields[5]:resid, fields[6-8]:X,Y,Z coordinates
            fields[3]=fields[3].strip() #if it is a standard aa with 3 letters, eliminate whitespace.
            fields[5]=fields[5].strip() #if it is a standard RESID with 4 characters, eliminate whitespace.
            #~ if fields[5]==pfields[5] and fields[3]!=pfields[3]: #avoids that same resid is used by different resnames.
                #~ return 'Corrupted PDB in position: '+pfields[5]+' Same resid has two or more different aminoacid codes/resnames'
            i=3
            while i<9:
                if fields[i].strip()=='':
                    return 'Missing required field in the PDB file at line: '+line
                i+=1

            if fields[5]!=pfields[5]: #resid has changed->new aa
                
                if fields[4]!=pfields[4]  or fields[9]!=pfields[9] or fields[5]=='1': #resid count has been reseted by new chain, new segid or whatever. 
                    ppos='0'
                    flag=0
                cpos=fields[5] #current position (resid) in the pdb during the present loop cycle
                if flag==1:
                    cpos2=int(str(cpos),16)
                    ppos2=int(str(ppos),16)
                elif flag==0:
                    cpos2=int(cpos)
                    ppos2=int(ppos)
                if cpos=='2710' and ppos=='9999':
                    cpos2=int(cpos,16)
                    flag=1
                    hexflag=1
                if (fields[4]==chain or chain == '') and cpos2 >= start and cpos2 <= stop and (segid in line[72:77] or segid==''):
                    if cpos2>=ppos2+1:
                        try:
                            seqplain.append([d[fields[3]],int(cpos)])
                        except: #Modified aminoacid
                            seqplain.append(('X',int(cpos)))

                    elif cpos2<ppos2 and cpos!=1: 
                        return 'Residue numbering order is corrupted in position:'+str(cpos2)

            pchain=fields[4]
            pfields=fields
            ppos=cpos

    fpdb.close()
    onlyaa='' #ignore the resids, pick the aa and that is it.
    for minilist in seqplain:
        onlyaa=onlyaa+minilist[0]
    #print(seqplain,onlyaa)
    if len(onlyaa)==0:
        return 'Unable to extract sequence from PDB file. Double check if the elements that define your interval exist: chain, segid, resid.'
    return (seqplain,onlyaa,hexflag)

def extract_PDB_seq_reading(pdb_filepath,model_chain_li):
    chain_seq_li=[]
    for chain_name in model_chain_li:
        thischain_data={}
        thischain_data["chainname"]=chain_name
        checkpdb_res=checkpdb_ngl(pdb_filepath, segid="",start=-1,stop=9999999999999999999, chain=chain_name)
        if not isinstance(checkpdb_res, tuple):
            raise Exception("Error extracting PDB sequence from the model")
        tablepdb,pdb_sequence,hexflag=checkpdb_res
        thischain_data["tablepdb"]=tablepdb
        thischain_data["pdbseq"]=pdb_sequence
        chain_seq_li.append(thischain_data)
    return chain_seq_li

def matchpdbfa_custom(uniprot_sequence,pdb_sequence, tablepdb):
    bestalig=pairwise2.align.localms(uniprot_sequence, pdb_sequence,100,-1,-10,-10)[0] #5,-1,-1.5,-1.5
    pdb_up_match=[]
    upalig=bestalig[0]
    pdbalig=bestalig[1]

    alig_pos=0
    up_resid=1 #we start at resid 1
    pdbseq_i=0
    while alig_pos < len(upalig):
        #[[ "PDB aa", "PDB resid"] , ["Up Seq AA", "UP Seq id"]]
        
        if pdbalig[alig_pos]=='-':
            aligpair=[["-",None],[upalig[alig_pos],up_resid]]
            pdb_up_match.append(aligpair)
            pdbseq_i-=1
        elif upalig[alig_pos]=='-':
            aligpair=[tablepdb[pdbseq_i],[upalig[alig_pos],None]]
            pdb_up_match.append(aligpair)            
            up_resid-=1
        else:
            aligpair=[tablepdb[pdbseq_i],[upalig[alig_pos],up_resid]]
            pdb_up_match.append(aligpair)
        pdbseq_i+=1
        up_resid+=1
        alig_pos+=1
    return pdb_up_match


def pdb_chain_to_uniprot_id(pdbid):
    fix_missing={'6LU7':['P0DTD1'],
     '6WPS': ['P0DTC2'],
     '6WPT': ['P0DTC2'],
     '6XCM': ['P0DTC2'],
     '6XCN': ['P0DTC2'],
     '6ZDH': ['P0DTC2'],
     '7K8T': ['P0DTC2'],
     '7K8U': ['P0DTC2'],
     '7K8W': ['P0DTC2'],
     '7K8X': ['P0DTC2'],
     '7K8Y': ['P0DTC2']}

    pdb_url="https://data.rcsb.org/rest/v1/core/entry/%s"%pdbid
    pdb_data=requests.get(pdb_url).json()
    pdb_polymers=pdb_data["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    chain_to_uniprot=[]
    for entityid in pdb_polymers: 
        entity_url="https://data.rcsb.org/rest/v1/core/polymer_entity/%s/%s"%(pdbid,entityid)
        entity_data=requests.get(entity_url).json()
        pdb_chains=entity_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
        #entity_data["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"]
        if "reference_sequence_identifiers" not in entity_data["rcsb_polymer_entity_container_identifiers"] and pdbid in fix_missing:
            uniprot_ids=fix_missing[pdbid]
        else:    
            uniprot_ids=[e["database_accession"] for e in entity_data["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"]if e['database_name']== 'UniProt']
        chain_to_uniprot.append({"uniprot_ids":uniprot_ids,"pdb_chains":pdb_chains})
    return chain_to_uniprot

def get_chain_model_match(best_chain_align,uniprot_chains,model_tablepdb):
    modelalig=best_chain_align["results"][0]
    chainalig=best_chain_align["results"][1]
    chainname=best_chain_align["chainname"]
    rcsb_chainseq=uniprot_chains[chainname]["chainseq"]
    chain_tablepdb=uniprot_chains[chainname]["tablepdb"]

    chain_model_match=[]
    alig_pos=0
    chain_seq_i=0
    model_seq_i=0
    while alig_pos < len(modelalig):
        #[[ "PDB chain aa", "PDB chain resid"] , ["Up Seq AA", "UP Seq id"]]
        
        if chainalig[alig_pos]=='-':
            aligpair=[["-",None],model_tablepdb[model_seq_i]]
            chain_model_match.append(aligpair)
            chain_seq_i-=1
        elif modelalig[alig_pos]=='-':
            aligpair=[chain_tablepdb[chain_seq_i],["-",None]]
            chain_model_match.append(aligpair)            
            model_seq_i-=1
        else:
            aligpair=[chain_tablepdb[chain_seq_i],model_tablepdb[model_seq_i]]
            chain_model_match.append(aligpair)
        chain_seq_i+=1
        model_seq_i+=1
        alig_pos+=1
    return chain_model_match

def check_for_missmatch_protportion(chain_model_match):
    totalcount=0
    mm_count=0
    for (chainpos,modelpos) in chain_model_match:
        if not chainpos[1]: #Not in the PDB chain related to the uniprot id
            continue
        if not modelpos[1]: #insertion in the model - doesn't relate to the uniprot seq
            continue
        totalcount+=1
        if chainpos[0] != modelpos[0]:
            mm_count+=1
    mm_prop=(mm_count/totalcount)*100
    return mm_prop

def check_avg_segment_len(chain_model_match):
    """
    Checks the average lengths of the segments between gaps- Takes the minimin lengths betweem the two alignmed sequences
    """
    n_segments_rcsb=0
    n_segments_model=0

    sum_seglen_rcsb=0
    sum_seglen_model=0

    len_last_segm_rcsb=0
    len_last_segm_model=0
    for (chainpos,modelpos) in chain_model_match:
        if chainpos[1]:
            len_last_segm_rcsb+=1
        else:
            if len_last_segm_rcsb:
                n_segments_rcsb+=1
                sum_seglen_rcsb+=len_last_segm_rcsb
                len_last_segm_rcsb=0
        if modelpos[1]:
            len_last_segm_model+=1
        else:
            if len_last_segm_model:
                n_segments_model+=1
                sum_seglen_model+=len_last_segm_model
                len_last_segm_model=0
    n_segments_rcsb+=1
    sum_seglen_rcsb+=len_last_segm_rcsb
    n_segments_model+=1
    sum_seglen_model+=len_last_segm_model
    
    seg_l_avg_rcsb=sum_seglen_rcsb/n_segments_rcsb
    seg_l_avg_model=sum_seglen_model/n_segments_model
    return min([seg_l_avg_rcsb,seg_l_avg_model])




def check_for_bad_alignment(chain_model_match,model_chainname):
    mm_prop=check_for_missmatch_protportion(chain_model_match)
    avg_seg_len=check_avg_segment_len(chain_model_match)
    if mm_prop > 30:
        print("Warning - Alignment rscb vs model with very high missmatches")
        if avg_seg_len < 3:
            print("Warning - Alignment rscb vs model with very short average lengths of segments between gaps")
            return False
    else:
        if avg_seg_len < 3:
            print("Warning - Alignment rscb vs model with very short average lengths of segments between gaps")
    return True
        


#basepath="/home/mariona/Documents/PhD/Research_projects/covid19/covid_platform_genomes/Data"
basepath="/protwis/sites/files/Covid19Data/Data"
tmp_pdb_path="/protwis/sites/files/Covid19Dynamics/tmp"
prot_data_path=os.path.join(basepath,"all_prot_data_filtered2_uniprot.data")

with open(prot_data_path, 'rb') as filehandle:  
    prot_data_up = pickle.load(filehandle)


protres_path=os.path.join(basepath,"protres.data")
with open(protres_path, 'rb') as filehandle:  
    d = pickle.load(filehandle)


with open("/gpcr/users/mariona/input/covid_dyn_dict.data", 'rb') as filehandle:  
    dyn_data = pickle.load(filehandle)




alldyndata=[]
for dyn_id, dynobj in dyn_data.items():
    try:
        pdb_filepath=dynobj["files"]["pdb"][0]["path"]
        pdbid=dynobj["pdbid"]
    except:
        print("COVID Dyn ID %s - no data" % dyn_id)
        continue
    if not os.path.isfile(pdb_filepath):
        print("PDB file not found: %s" % pdb_filepath)
        continue
    thisdyndata={}
    thisdyndata["dyn_id"]=dyn_id
    thisdyndata["pdbid"]=pdbid
    thisdyndata["pdb_filepath"]=pdb_filepath
    alldyndata.append(thisdyndata)

del dyn_data
gc.collect()

overwrite=False

final_dyn_data={}
out_path=os.path.join(basepath,"mapping_results.data")
if (not overwrite) and os.path.isfile(out_path):
    with open(out_path, 'rb') as filehandle:  
        final_dyn_data = pickle.load(filehandle)



map_pdb_to_up={}
for dynobj in alldyndata:
    dyn_id=dynobj["dyn_id"]
    if dyn_id in final_dyn_data:
        continue
    print("Dyn ID: %s" % dyn_id)

    try:
        pdbid=dynobj["pdbid"]
        pdb_filepath=dynobj["pdb_filepath"]
        url="https://data.rcsb.org/rest/v1/core/entry/%s"%pdbid
        pdb_data=requests.get(url).json()

        #1) Map variants to RCSB PDB file

        # I need a dict relating Uniprot ID in the PDB to chain in the PDB
        try:
            chain_to_uniprot=pdb_chain_to_uniprot_id(pdbid)
        except:
            print("No PDB data for PDB ID %s"%pdbid)
            continue

        chain_to_uniprot_ordered={}
        for c_u_set in chain_to_uniprot:
            myuniprot_ids=c_u_set["uniprot_ids"]
            mypdb_chains=c_u_set["pdb_chains"]
            myuniprot_ids_tp=tuple(sorted(myuniprot_ids))
            if myuniprot_ids_tp not in chain_to_uniprot_ordered:
                chain_to_uniprot_ordered[myuniprot_ids_tp]=set()
            current=chain_to_uniprot_ordered[myuniprot_ids_tp]
            mypdb_chains_s=set(mypdb_chains)
            chain_to_uniprot_ordered[myuniprot_ids_tp]=current.union(mypdb_chains_s)

        #Download RCSB PDB file
        rawpdb_path=os.path.join(tmp_pdb_path,"%s.pdb"%pdbid)

        if not os.path.isfile(rawpdb_path):
            urllib.request.urlretrieve('https://files.rcsb.org/download/%s.pdb'%pdbid, rawpdb_path)
        #Extract sequence from each chain of the RCSB PDB file
        raw_chain_seq_li=extract_PDB_seq(rawpdb_path)
        raw_chain_seq={e["chainname"]:e for e in raw_chain_seq_li}

        final_dyn_data[dyn_id]={}
        final_dyn_data[dyn_id]["pdbid"]=pdbid
        final_dyn_data[dyn_id]["pdb_filepath"]=pdb_filepath

        for (uniprot_ids,pdb_chains) in chain_to_uniprot_ordered.items(): #iterate relation uniprot id - pdb chains 
            uniprot_ids_ok=[up for up in uniprot_ids if (up in prot_data_up)]
            if uniprot_ids_ok:
                if len(uniprot_ids_ok)>1:
                    raise Exception("Dyn %s - More than one uniprot ID for PDB %s chainset %s : %s"%(dyn_id,pdbid,",".join(pdb_chains),",".join(uniprot_ids_ok)))
                #Uniprot ID of Sars-Cov-2 to be aligned to the corresponding chains in RCSB PDB file
                uniprot_id=uniprot_ids_ok[0]
                uniprot_sequence=prot_data_up[uniprot_id]["seq"]
                up_to_rcsbchains={(p,i+1):[] for i,p in enumerate(uniprot_sequence)}
                up_to_model=copy.deepcopy(up_to_rcsbchains)
                map_pdb_to_up[pdbid]={}
                uniprot_chains={}
                for chain in pdb_chains:
                    map_pdb_to_up[pdbid][chain]={}
                    map_pdb_to_up[pdbid][chain]["uniprotid"]=uniprot_id
                    chaininfo=raw_chain_seq[chain]
                    tablepdb=chaininfo["tablepdb"]
                    pdb_sequence=chaininfo["pdbseq"]
                    pdb_up_match=matchpdbfa_custom(uniprot_sequence,pdb_sequence, tablepdb)
                    rcsbpdb_to_up={}
                    for (rcsbpdbpos,uppos) in pdb_up_match:
                        if rcsbpdbpos[1]:
                            rcsbpdb_to_up[tuple(rcsbpdbpos)]=uppos
                        if uppos[1]:
                            uppost=tuple(uppos)
                            if uppost in up_to_rcsbchains:
                                up_to_rcsbchains[uppost].append(chain)
                    uniprot_chains[chain]={"chainseq":pdb_sequence,"rcsbpdb_to_up":rcsbpdb_to_up,"tablepdb":tablepdb}
                    map_pdb_to_up[pdbid][chain]["rcsbpdb_to_up"]=rcsbpdb_to_up
                #check if the different chains corresponding to the uniprot ID are equivalent or cover different parts of the uniprot seq. Usually in a RCSB PDB file we have 1 chain or several equivalent chains for each uniprot ID 
                pos_corresponding=len([e for e in up_to_rcsbchains.values() if len(e)==len(pdb_chains)])
                ccorresp_score=pos_corresponding/len(uniprot_sequence)
                if ccorresp_score < 0.5:
                    raise Exception("RCSB PDB chains corresponding to uniprot ID %s are not redundant"% uniprot_id)
                #2) Align RCSB PDB file to model
                try:
                    model_chain_seq_li=extract_PDB_seq(pdb_filepath,True)
                except IndexError:
                    print("Dyn %s - Alternative!!!" % dyn_id)
                    model_chain_li=obtain_prot_chains(pdb_filepath)
                    model_chain_seq_li = extract_PDB_seq_reading(pdb_filepath,model_chain_li)

                #The model PDB chain may correspond to 1 RCSB PDB chain (and thus uniprot ids) or multiple. For each model chain I need to align 1 or more RCSB PDB chain
                for model_chaininfo in model_chain_seq_li: #all model pos
                    model_tablepdb=model_chaininfo["tablepdb"]
                    model_pdb_sequence=model_chaininfo["pdbseq"]
                    model_chainname=model_chaininfo["chainname"]
                    best_chain_align=None
                    for chainname, rcsb_chaininfo in uniprot_chains.items():
                        #Since we may have several RCSB PDB chain for this uniprot ID (presumably equivalent), we will keep the one with higher score
                        rcsb_chainseq=rcsb_chaininfo["chainseq"]
                    
                        bestalig=pairwise2.align.localms(model_pdb_sequence, rcsb_chainseq,100,-1,-10,-10)[0] #5,-1,-1.5,-1.5
                        if not best_chain_align or bestalig[2]>best_chain_align["results"][2]:
                            best_chain_align={"results":bestalig[0:3],"chainname":chainname}

                    chain_model_match=get_chain_model_match(best_chain_align,uniprot_chains,model_tablepdb)
                    if len(model_chain_seq_li) >1:
                        keep_chain=check_for_bad_alignment(chain_model_match,model_chainname) #Can be than a model PDB chain corresponds to a RCSB PDB chain that has been discarded for not being part of a covid protein. If that happens, the alignment will be very bad. 
                        if not keep_chain:
                            print("Warning - Ignoring model chain %s." % model_chainname)
                            continue
                    rcsbpdb_to_up=uniprot_chains[best_chain_align["chainname"]]["rcsbpdb_to_up"]
                    #up_to_model
                    for (chainpos,modelpos) in chain_model_match:
                        if not chainpos[1]: #Not in the PDB chain related to the uniprot id
                            continue
                        if not modelpos[1]: #insertion in the model - doesn't relate to the uniprot seq
                            continue
                        up_pos=rcsbpdb_to_up[tuple(chainpos)]
                        ca_index=None
                        if len(modelpos)>2:
                            ca_index=modelpos[2]
                        if tuple(up_pos) in up_to_model:
                            up_to_model[tuple(up_pos)]={"aa":modelpos[0],"pos":modelpos[1],"chain":model_chainname, "ca_index":ca_index}
                if "up_to_model" in final_dyn_data[dyn_id]:
                    raise Exception("More than one uniprot-model mapping found in dyn %s" % dyn_id)
                final_dyn_data[dyn_id]["up_to_model"]=up_to_model
                final_dyn_data[dyn_id]["uniprot_id"]=uniprot_id
    except Exception as e:
        print("ERROR in dyn %s: %s" %(dyn_id,e))        

out_path=os.path.join(basepath,"mapping_results.data")
print("Save file: ",out_path)
with open(out_path, "wb") as filehandle:  
    # store the data as binary data stream
    pickle.dump(final_dyn_data, filehandle)


out_path2=os.path.join(basepath,"mapping_results_rcsb_up.data")
print("Save extra file: ",out_path2)
with open(out_path2, "wb") as filehandle:  
    # store the data as binary data stream
    pickle.dump(map_pdb_to_up, filehandle)

#map_path=os.path.join(basepath,"mapping_results.data")
#with open(map_path, "rb") as filehandle:  
#    # store the data as binary data stream
#    final_dyn_data=pickle.load(filehandle)

