# Run in vagrant/server

import os
import pickle
import pandas as pd
import re
from covid19.models import *
from django.core.exceptions import MultipleObjectsReturned
import requests
from io import StringIO

def query_uprot_singleentry(uniprotkbac,fields):
    payload = {'query': 'id:%s'%uniprotkbac,'format': 'tab','columns': ",".join(fields)}
    result = requests.get("http://www.uniprot.org/uniprot/", params=payload)
    if result.ok and result.text:
        result_str=result.text
        (headers,vals)=result_str.split("\n")[:2]
        results_dict=dict(zip(headers.split("\t"),vals.split("\t")))
    else:
        results_dict={}
    return results_dict


basepath="/protwis/sites/files/Covid19Data/Data"
#1)Save domain info to DB
save_domains_to_db=True
if save_domains_to_db:
    domains_info_path=os.path.join(basepath,"domains.xlsx")
    domains_info=pd.read_excel(domains_info_path,index_col=0).fillna(value="")

    ok_columns=['Active site', 'Binding site', 'Calcium binding',
            'Cofactor', 'DNA binding',  'Kinetics', 'Metal binding', 
           'Nucleotide binding', 'Pathway', 'pH dependence', 'Redox potential',
           'Site', 'Temperature dependence', 'Coiled coil', 'Compositional bias',
            'Domain [FT]', 'Motif',  'Region',
           'Repeat', 'Zinc finger']
    fix_domnames={
        "Site": "Site of interest",
        "Domain [FT]": "Domain of interest",
        'Coiled coil':"coiled-coil region",
        'Metal binding':'Metal ion-binding site',
        "Nucleotide binding": "Nucleotide-bnding region",
        "Region":"Region of interest",
        'Zinc finger':'Zinc finger region'
    }
    alldata_save=[]
    for columnName in ok_columns:
        print('Colunm Name : ', columnName)
        columnData=domains_info[columnName]
        for uniprotid, prot_col_data in columnData.iteritems():
            up_seq=domains_info["Sequence"][uniprotid]
            if prot_col_data:
                regions_li=re.split(r';\s(?=(?:[^"]*"[^"]*")*[^"]*$)',prot_col_data)
                for region_str in regions_li:
                    region_data_all=re.split(r'\s/(?=(?:[^"]*"[^"]*")*[^"]*$)',region_str)
                    isfirst=True
                    region_save={}
                    for region_data in region_data_all:
                        if isfirst:
                            pos_data=region_data.split(" ")[1]
                            if ".." in pos_data:
                                (frompos,topos)=[int(e) for e in pos_data.split("..")]
                            else:
                                frompos=int(pos_data)
                                topos=frompos
                            region_name=columnName
                            if region_name in fix_domnames:
                                region_name=fix_domnames[region_name]
                            region_save["region_name"]=region_name
                            region_save["frompos"]=frompos
                            region_save["topos"]=topos
                            region_save["uniprotkbac"]=uniprotid
                            region_save["dom_seq"]=up_seq[frompos-1:topos]
                            isfirst=False
                        else:
                            (fieldname,fieldval)=re.split(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)',region_data)
                            region_save[fieldname] = fieldval.replace('"','')
                    alldata_save.append( region_save)
    for tosave in alldata_save:
        domnote=tosave.get("note","")
        if domnote:
            if "; structural" in domnote:
                (dom1,dom2)=domnote.split("; ")
                domnote="%s (%s)"%(dom1,dom2)
        domobj, created_domobj = CovidDomains.objects.get_or_create(
            frompos = tosave["frompos"],
            topos = tosave["topos"],
            region_name = tosave["region_name"],
            uniprotkbac = tosave["uniprotkbac"],
            note = domnote,
            evidence = tosave.get("evidence",""),
            domain_seq=tosave["dom_seq"]
        )
else:
    print("save_domains_to_db is set to False - no domains are saved")


#2) Save mapping info to DB
save_mapping_to_db=True
if save_mapping_to_db:
    mapping_info_path=os.path.join(basepath,"mapping_results.data")
    with open(mapping_info_path, 'rb') as filehandle:  
        dyn_up_to_model = pickle.load(filehandle)


    for dyn_id,dyn_data in dyn_up_to_model.items():
        try:
            uniprotid=dyn_data["uniprot_id"]
            up_to_model=dyn_data["up_to_model"]
        except KeyError:
            print("Missing data for dyn %s"% dyn_id)
            continue            
        try:
            modelfileobj=CovidFiles.objects.get(id_file_types__is_model=True,covidfilesdynamics__id_dynamics=dyn_id)
        except:
            print("No model file found for dyn %s"% dyn_id)
            continue
        try:
            protobj = CovidProtein.objects.get(uniprotkbac=uniprotid)
        except (CovidProtein.DoesNotExist):
            fields=["entry_name","protein_names","organism"]
            results_dict=query_uprot_singleentry(uniprotid,fields)
            if results_dict:
                print("Saving new CovidProtein for %s"% uniprotid)
                prot_names_all=results_dict["Protein names"]
                prot_name=prot_names_all[:prot_names_all.find("(")].strip()
                uniprot_entry = results_dict["Entry name"]
                species =results_dict["Organism"]
                protobj, created_protobj = CovidProtein.objects.get_or_create(
                    uniprotkbac=uniprotid,
                    uniprot_entry=uniprot_entry,
                    name=prot_name,
                    species=species
                )
            else:
                print("uniprot ID not found")
                continue
        except MultipleObjectsReturned:
            print("More than one protein entry for %s (dyn %s)"%(uniprotid,dyn_id))
            continue
        for up_pos, model_pos in sorted(up_to_model.items(),key=lambda x: x[0][1]):
            #('L', 461)
            upposobj, created_upposobj = CovidUniprotSeqPositions.objects.get_or_create(
                    id_protein=protobj,
                    seqpos=up_pos[1],
                    aa=up_pos[0]
                )
            if model_pos:
                #{'aa': 'L', 'pos': 359, 'chain': 'B'}
                ca_index=model_pos.get("ca_index")
                modelposobj, created_modelposobj = CovidModelSeqPositions.objects.get_or_create(
                        id_file = modelfileobj,
                        seqpos =model_pos["pos"],
                        aa =model_pos["aa"],
                        chainid = model_pos["chain"],
                        id_uniprotpos=upposobj,
                        ca_atom_index=ca_index
                    )



#3) Get seq positions of all final COVID-19 proteins relative to the uniprot sequence

protname_map={
 "2'-O-methyltransferase": None,
 '3C-like proteinase':"NSP5" ,
 'Envelope small membrane protein':"E" ,
 'Helicase': "NSP13",
 'Host translation inhibitor nsp1':"NSP1" ,
 'Membrane protein': "M",
 'Non-structural protein 10': "NSP10",
 'Non-structural protein 11': "NSP11",
 'Non-structural protein 2': "NSP2",
 'Non-structural protein 3': "NSP3",
 'Non-structural protein 4': "NSP4",
 'Non-structural protein 6': "NSP6",
 'Non-structural protein 7': "NSP7",
 'Non-structural protein 8': "NSP8",
 'Non-structural protein 9': "NSP9",
 'Nucleoprotein': "N",
 'ORF3a protein': "ORF3a",
 'ORF6 protein': "ORF6",
 'ORF7a protein': "ORF7a",
 'ORF7b protein': "ORF7b",
 'ORF8 protein': "ORF8",
 'ORF9b protein': "ORF9b",
 'Proofreading exoribonuclease': "NSP14",
 'RNA-directed RNA polymerase': "NSP12",
 'Replicase polyprotein 1a': None,
 'Replicase polyprotein 1ab': None,
 'Spike glycoprotein': "Spike",
 'Spike protein S1': None,
 'Spike protein S2': None,
 "Spike protein S2'": None,
 'Uncharacterized protein 14':"ORF14" ,
 'Uridylate-specific endoribonuclease':"NSP15"
 }


up_to_protpos={}
for protobj in CovidProtein.objects.filter(species__icontains="SARS-CoV-2"):
    uniprot_id=protobj.uniprotkbac
    query_res=query_uprot_singleentry(uniprot_id,["feature(CHAIN)","sequence"])
    if uniprot_id in up_to_protpos:
        continue
    
    if query_res:
        prot_ptm_chains=query_res["Chain"]
        up_seq=query_res["Sequence"]
        ptm_chain_li=prot_ptm_chains.split("; CHAIN ")
        for ptm_chain_s in ptm_chain_li:
            ptm_chain_s=ptm_chain_s.lstrip("CHAIN ")
            ptm_chain_allinfo=re.split(r';\s(?=(?:[^"]*"[^"]*")*[^"]*$)',ptm_chain_s)
            isfirst=True
            finalprot_data={}
            save_finalprot_data=False
            if ptm_chain_allinfo==[""]:
                continue
            for ptm_chain_info in ptm_chain_allinfo:
                if isfirst:
                    isfirst=False
                    if ".." in ptm_chain_info:
                        (frompos,topos)=[int(e) for e in ptm_chain_info.split("..")]
                    else:
                        frompos=int(ptm_chain_info)
                        topos=frompos
                    finalprot_data["from"]=frompos
                    finalprot_data["to"]=topos
                else:
                    (fieldname,fieldval)=re.split(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)',ptm_chain_info)
                    fieldname=fieldname.lstrip(" /")
                    fieldval=fieldval.strip('"')
                    if fieldname=="note":
                        if fieldval in protname_map and protname_map[fieldval]:
                            finalprot_data["prot"]=protname_map[fieldval]
                            save_finalprot_data=True
                    #if fieldname not in ["id","note"]:
                    #    print(fieldname,fieldval)
            if save_finalprot_data:
                if uniprot_id not in up_to_protpos:
                    up_to_protpos[uniprot_id]=[]
                up_to_protpos[uniprot_id].append(finalprot_data)
                if finalprot_data["prot"]=="Spike":
                    finalprot_data["from"]=1 #For some reason uniprot says the final protein starts at pos 13!
                finprot_obj,created_finprot_obj=CovidFinalProtein.objects.get_or_create(name=finalprot_data["prot"])
                protfinprot_obj,created_protfinprot_obj=CovidProteinFinalprotein.objects.get_or_create(
                    id_protein=protobj,
                    id_finalprotein=finprot_obj
                    )
                protfinprot_obj.finalprot_seq_start=finalprot_data["from"]
                protfinprot_obj.finalprot_seq_end=finalprot_data["to"]
                protfinprot_obj.save()
                #print(">"+finalprot_data["prot"])
                #print(up_seq[finalprot_data["from"]-1:finalprot_data["to"]])





else:
    print("save_domains_to_db is set to False - mapping not saved")