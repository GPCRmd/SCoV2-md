
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import pickle
import requests
from io import StringIO
import pandas as pd

def query_uprot(query,fields):
    payload = {'query': query,'format': 'tab','columns': ",".join(fields)}
    result = requests.get("http://www.uniprot.org/uniprot/", params=payload)
    if result.ok and result.text:
        result_str=result.text
        upresults=pd.read_csv(StringIO(result_str), sep='\t', header=0, index_col=0)
        return upresults
    else:
        return False


def classify_seq_data(basepath,print_results=False):
    #Protein sequences downloaded from NCBI 
    #https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049
    orig_protseq=os.path.join(basepath,"all_protein_sequences.fasta")

    #Uniprot data
    up_data_path=os.path.join(basepath,"sarscov2_uniprot.csv")
    if os.path.isfile(up_data_path):
        up_data=pd.read_csv(up_data_path, header=0, index_col=0)
    else:
        myquery="organism:2697049"
        up_data=query_uprot(myquery,['id','entry_name', 'protein_names', 'reviewed', 'keywords', 'database(RefSeq)'])
        up_data.to_csv(up_data_path)

    prot_data={}
    #Group seq. data by protein, differentiating ref seq. to other sequenced genomes
    for seq_record in SeqIO.parse(orig_protseq, 'fasta'):
        is_refSeq=False
        if seq_record.id.startswith("YP_"):
            is_refSeq=True
        description=seq_record.description
        description=description.split("|")[1]
        prot_name=description.split("[")[0]
        prot_name=prot_name.strip()
        prot_name=prot_name.replace(", partial","")
        if prot_name=="ORF7b": #some entries are ORF7b and other ORF7b protein
            prot_name="ORF7b protein"
        if prot_name not in prot_data:
            prot_data[prot_name]={}
            prot_data[prot_name]["var_seqs"]=[]
            prot_data[prot_name]["uniprot_data"]={}
        read_d={}
        read_d["id"]=seq_record.id
        read_d["description"]=seq_record.description
        read_d["seq"]=seq_record.seq
        if is_refSeq:
            if "refSeq" in prot_data[prot_name]:
                prot_data[prot_name]["refSeq"]["id"].append(read_d["id"])
            else:
                read_d["id"]=[read_d["id"]]
                prot_data[prot_name]["refSeq"]=read_d
        else:
            prot_data[prot_name]["var_seqs"].append(read_d)

    #Assign uniprot ID to each protein based on the protein name
    for prot_name,myprot in prot_data.items():
        if myprot["var_seqs"]:
            refid=myprot["refSeq"]["id"][0]+";"
            myupdata=up_data.loc[(up_data['Cross-reference (RefSeq)']==refid) & (up_data["Status"]=="reviewed") ]
            if myupdata.empty:
                myupdata=up_data.loc[(up_data['Cross-reference (RefSeq)']==refid) & (up_data["Keywords"]=="Reference proteome")]
            if len(myupdata)==1:
                uprotid=myupdata.index[0]
                myprot["uniprot_data"]["uniprotid"]=uprotid
            elif len(myupdata)>1:
                raise Exception("More than one match for protein %s"% prot_name)

    if print_results:
        print("Prot name\tTotal seqs.\tnon-refSeq")          
        for k in sorted(prot_data.keys()):
            v=prot_data[k]                                   
            print(k,"\t",len(v["var_seqs"])+len(v["refSeq"]["id"]),"\t",len(v["var_seqs"]))
    return prot_data



def align_ref_var(refseq,varseq):
    bestalig=pairwise2.align.localms(refseq,varseq,5,-1,-1.5,-1.5)[0]
    return bestalig

def obtain_vars(refseq,varseq,this_varprot,thisprot):
    bestalig=align_ref_var(refseq,varseq)
    (aligned_refseq,aligned_varseq)=bestalig[:2]
    i=0
    seqpos=0
    protvars=[]
    for refpos in aligned_refseq:
        varpos=aligned_varseq[i]
        if refpos!=varpos:
            if refpos=="-":
                if seqpos==0:
                    raise Exception("Deletion of 1st residue not considered")
                seqpos-=1
                refpos=refseq[i-1]
                if refpos=="-":
                    raise Exception("Deletion of more than one position not considered")
                varpos=refpos+varpos

            var_info={"from":refpos,"to":varpos,"position":seqpos}
            protvars.append(var_info)
            varname="%s;%s;%s"%(seqpos,refpos,varpos)
            if varname not in thisprot["all_found_variants"]:
                thisprot["all_found_variants"][varname]=[]
            thisprot["all_found_variants"][varname].append(this_varprot["id"])
        elif refpos=="-": #both are "-"
            seqpos-=1
        seqpos+=1
        i+=1
    return protvars

def obtain_all_prot_variants(prot_data):
    for protname in sorted(prot_data.keys()):
        thisprot=prot_data[protname]
        if thisprot["var_seqs"]:
            refseq=thisprot['refSeq']['seq']
            varprots=thisprot["var_seqs"]
            thisprot["all_found_variants"]={}
            for this_varprot in varprots:
                varseq=this_varprot["seq"]
                if varseq!=refseq:
                    protvars=obtain_vars(refseq,varseq,this_varprot,thisprot)
                    if not type(protvars)==list:
                        break
                    this_varprot["found_variants"]=protvars

    return prot_data

def apply_filter(aato):
    if aato!="X":
        return True 
    else:
        return False
        
def remove_from_general_prot_vars(thisprot,varprot_var_info,var_id):
    if "all_found_variants" in thisprot:
        allvars=thisprot["all_found_variants"] 
        varname="%s;%s;%s"%(varprot_var_info["position"],varprot_var_info["from"],varprot_var_info["to"])
        if varname in allvars:
            allvars[varname]
            if var_id in allvars[varname]:
                allvars[varname].remove(var_id)
                if len(allvars[varname])==0:
                    del allvars[varname]
                if len(allvars)==0:
                    del thisprot["all_found_variants"] 

def filter_out_partial_sequence_effects(prot_data):
    #remove "X" residues
    for protname in sorted(prot_data.keys()):
        thisprot=prot_data[protname]
        if "all_found_variants" in thisprot:
            allvars=thisprot["all_found_variants"] 
            new_allvars={}
            for varname,seqids in allvars.items():
                (pos,aafrom,aato)=varname.split(";")
                passed_filter=apply_filter(aato)
                if passed_filter:
                    new_allvars[varname]=seqids
            thisprot["all_found_variants"]=new_allvars

            varprots=thisprot["var_seqs"]
            for this_varprot in varprots:
                if "found_variants" in this_varprot:
                    varprot_vars=this_varprot["found_variants"]
                    new_varprot_vars=[]
                    for varprot_var_info in varprot_vars:
                        aato=varprot_var_info["to"]
                        passed_filter=apply_filter(aato)
                        if passed_filter:
                            new_varprot_vars.append(varprot_var_info)
                    this_varprot["found_variants"]=new_varprot_vars

    #Remove false deletions caused by partial sequence
    for protname in sorted(prot_data.keys()):
        thisprot=prot_data[protname]
        if thisprot["var_seqs"]:
            refseq=thisprot['refSeq']['seq']
            varprots=thisprot["var_seqs"]
            for this_varprot in varprots:
                description=this_varprot["description"]
                var_id=this_varprot["id"]
                is_partial="partial" in description
                if "found_variants" in this_varprot:
                    if is_partial:
                        varprot_vars=sorted(this_varprot["found_variants"], key=lambda x:x["position"])
                        if varprot_vars[0]["position"]==0 and varprot_vars[0]["to"]=="-":
                            new_varprot_vars=[]
                            skip_var=True
                            for i,varprot_var_info in enumerate(varprot_vars):
                                if varprot_var_info["position"]!=i or varprot_var_info["to"]!="-":
                                    skip_var=False
                                if skip_var:
                                    remove_from_general_prot_vars(thisprot,varprot_var_info,var_id)
                                else:
                                    new_varprot_vars.append(varprot_var_info)
                            if new_varprot_vars:
                                this_varprot["found_variants"]=new_varprot_vars
                            else:
                                del this_varprot['found_variants']
                        last_pos=len(refseq)-1
                        if varprot_vars[-1]["position"]==last_pos and varprot_vars[-1]["to"]=="-":
                            new_varprot_vars=[]
                            skip_var=True
                            for i,varprot_var_info in enumerate(varprot_vars[::-1]):
                                negi=last_pos-i
                                if varprot_var_info["position"]!=negi or varprot_var_info["to"]!="-":
                                    skip_var=False
                                if skip_var:
                                    remove_from_general_prot_vars(thisprot,varprot_var_info,var_id)
                                else:
                                    new_varprot_vars.append(varprot_var_info)
                            if new_varprot_vars:
                                this_varprot["found_variants"]=new_varprot_vars
                            else:
                                del this_varprot['found_variants']
    return prot_data


basepath="/home/mariona/Documents/PhD/Research_projects/covid19/covid_platform_genomes/Data"
#Load data from NCBI and classify both ref. seq and other sequences by protein
prot_data=classify_seq_data(basepath)

#Align variant genomes to refSeqs to obtain variants
prot_data=obtain_all_prot_variants(prot_data)
out_dict_path=os.path.join(basepath,"all_prot_data.data")
with open(out_dict_path, "wb") as filehandle:  
    # store the data as binary data stream
    pickle.dump(prot_data, filehandle)

# Filter variant data to not count as deletion the cases in which part of the protein is not sequenced 
prot_data=filter_out_partial_sequence_effects(prot_data)



out_dict_path_filt=os.path.join(basepath,"all_prot_data_filtered2.data")
with open(out_dict_path_filt, "wb") as filehandle:  
    # store the data as binary data stream
    pickle.dump(prot_data, filehandle)


#Change of format (in case SeqRecord doesn't work) - needed to run map_info_to_struc.py
prot_data_up={}
for protname,thisprot in prot_data.items():
    if "uniprot_data" in thisprot:
        myuniprotid=thisprot["uniprot_data"]["uniprotid"]        
        refseq=thisprot["refSeq"]
        refseq["seq"]=str(refseq["seq"])
        prot_data_up[myuniprotid]=refseq

out_dict_path_filt_up=os.path.join(basepath,"all_prot_data_filtered2_uniprot.data")
with open(out_dict_path_filt_up, "wb") as filehandle:  
    # store the data as binary data stream
    pickle.dump(prot_data_up, filehandle)