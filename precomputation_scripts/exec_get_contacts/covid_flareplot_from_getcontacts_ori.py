import argparse
import gc
import os
#import mdtraj as md
import json
import pickle
import requests
import re
from Bio.PDB.Polypeptide import *
import seaborn as sns

parser = argparse.ArgumentParser(description="""
Uses the oytput of the scrpts that generate the input for meta-analysis page (getcontacts) to generate json files for the flareplots in the Workbench.
""") 
parser.add_argument(
   '-i',
    dest='dyn_dict_path',
    default="/gpcr/users/mariona/input/covid_dyn_dict.data",
    action='store',
    type=str,
    help='Path to the input dictionary wit data on the simulations'
)
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
    '--ignore_publication',
    action='store_true',
    dest='ignore_publication',
    default=False,
    help='Consider both published and unpublished dynamics.',
)
parser.add_argument(
    '--overwrite',
    dest='overwrite',
    action='store_true',
    default=False,
    help='(bool) repeat get_dynamic_contacts even if there already exists a results file for this files'
)


def generate_data_dict():
    cont_li=['sb', 'pc', 'ps', 'ts', 'vdw', 'hb', 'wb', 'wb2', "hp"]
    traj_int_d={e:dict() for e in cont_li}
    for int_d in traj_int_d.values():
        int_d['defaults']={'edgeColor': 'rgba(50,50,50,100)', 'edgeWidth': 2}
        int_d["trees"]=[{
                "treeLabel" :'Helices',
                "treePaths":set(), # ex. '2.1x30'
            }]
        int_d["edges"]=dict() #list of {'name1':, 'frames':, 'name2':, 'helixpos':} -> helixpos can be Intra or Inter
                        #ex: {"frames": [ 0, 10, 12],  'helixpos': 'Intra', 'name1': '5x38', 'name2': '5x39'}
        int_d["tracks"]=[
                    {"trackLabel":"Helices",
                     "trackProperties" :list() # {'color': '#79C268', 'nodeName': '3x25', 'size': '1.0'}
                    }
                    ]
    return traj_int_d

def save_json(dyn_id,traj_id,int_type,int_data,resultspath):
    json_filename="%s_trj_%s_%s.json" % (traj_id,dyn_id,int_type)
    filpath=os.path.join(resultspath,int_type)
    if not os.path.isdir(filpath):
        os.makedirs(filpath)
    json_namepath=os.path.join(filpath,json_filename)
    print("Saved: %s" % json_namepath)
    with open(json_namepath, 'w') as outfile:
        json.dump(int_data, outfile)

def extract_res_info(res1,all_chains):
    res1_post=res1[:res1.rfind(":")] #remove atom info
    (chain,resname3,resid)=res1_post.split(":")
    resid=re.sub('\D', '', resid)
    try:
        resname=three_to_one(resname3)
    except:
        resname="X"
    gnum1= "%s%s_%s"%(resname,resid,chain) #name that will appear at flare plot
    if chain not in all_chains:
        all_chains.append(chain)
    h1=chain
    treep_num=all_chains.index(chain) +1
    treep1="%s.%s" % (treep_num,gnum1)
    return (gnum1,treep1,h1,int(resid))

def apply_colors(int_data):
    mytracks= int_data["tracks"][0]["trackProperties"] 
    min_val=False
    max_val=False
    first=True
    for e in mytracks:
        nodename=e["nodeName"]
        resid=int(nodename.split("_")[0][1:])
        if first:
            min_val=resid
            max_val=resid
            first=False
        else:
            if resid>max_val:
                max_val=resid
            elif resid < min_val:
                min_val=resid
    num_resids=max_val - min_val+1
    colormap=sns.color_palette("Spectral",n_colors=num_resids).as_hex()
    first=True
    for e in mytracks:
        nodename=e["nodeName"]
        resid=int(nodename.split("_")[0][1:])
        resid_pos=resid-min_val
        this_color=colormap[resid_pos]
        e["color"]=this_color
    return int_data


def create_p_jsons(dynfiles_traj,gcdata_path_dyn,dyn_id,resultspath):
    hb_list=["hbbb","hbsb","hbss"]
    for traj in dynfiles_traj:
        n_frames=traj["framenum"]
        traj_id=traj["id"]
        print("\tTraj id: %s"%traj_id)
        traj_int_d=generate_data_dict()

        gcdata_path_dyn_traj=os.path.join(gcdata_path_dyn,"dyn%(dynid)s-%(trajid)s_dynamic.tsv" % {"dynid" :dyn_id, "trajid":traj_id} )
        pre_frame=False
        all_chains=[]
        with open(gcdata_path_dyn_traj) as infile:
            for line in infile:
                line = line.strip()
                if "total_frames" in line:
                    el = line.split(" ")
                    file_total_frames=int(el[1].split(":")[1])
                if len(line) == 0 or line[0] == "#":
                    continue
                allinfo = line.split("\t")
                if len(allinfo)==4:
                    (frame,int_type,res1,res2)=allinfo
                elif len(allinfo)==5:
                    (frame,int_type,res1,res2,res3)=allinfo
                elif len(allinfo)==6:
                    (frame,int_type,res1,res2,res3,res4)=allinfo
                else:
                    print("Incorrect number of elements in line. Skipping. Line: %s"%line)
                    continue
#                if frame != pre_frame:
#                    if int(pre_frame) == frame_ends_bytraj[traj_rep]:#[!]
#                        print("\tTraj id: %s"%traj_id)
                #add all res:
                resinfo1=extract_res_info(res1,all_chains)
                resinfo2=extract_res_info(res2,all_chains)
                if resinfo1 and resinfo2:
                    (gnum1,treep1,h1,resid1)=resinfo1
                    (gnum2,treep2,h2,resid2)=resinfo2
                else:
                    continue
                for int_typeid, int_data in traj_int_d.items():
                    if treep1 not in int_data["trees"][0]["treePaths"]:
                        int_data["trees"][0]["treePaths"].add(treep1) # ex. '2.1x30'
                        #int_data["tracks"][0]["trackProperties"].append({'color': nodecolor1, 'nodeName': gnum1, 'size': '1.0'})
                        int_data["tracks"][0]["trackProperties"].append({'nodeName': gnum1, 'size': '1.0'})
                    if treep2 not in int_data["trees"][0]["treePaths"]:
                        int_data["trees"][0]["treePaths"].add(treep2)
                        #int_data["tracks"][0]["trackProperties"].append({'color': nodecolor2, 'nodeName': gnum2, 'size': '1.0'})
                        int_data["tracks"][0]["trackProperties"].append({ 'nodeName': gnum2, 'size': '1.0'})
                #add this particular inteaction
                if int_type in hb_list:
                    int_type="hb"
                if int_type in traj_int_d:
                    edge_d=traj_int_d[int_type]["edges"]
                    if (gnum1,gnum2) in edge_d:
                        edge_d[(gnum1,gnum2)]["frames"].append(frame)
                    elif (gnum2,gnum1) in edge_d:
                        edge_d[(gnum2,gnum1)]["frames"].append(frame)
                    else:
                        if (h1==h2):
                            hpos="Intra"
                        else:
                            hpos="Inter"
                        is_neighbour=False
                        if abs(resid1-resid2) <=5:
                            is_neighbour=True
                        edge_d[(gnum1,gnum2)] = {'name1':gnum1 , 'name2':gnum2 , 'frames':[frame],  'helixpos':hpos, "is_neighbour":is_neighbour}

                pre_frame=frame

        for int_type, int_data in traj_int_d.items():
            if not int_data["edges"]:
                print("No data for %s" % int_type)
                continue
            int_data["trees"][0]["treePaths"] = list(int_data["trees"][0]["treePaths"])
            int_data["edges"] = [v for k,v in int_data["edges"].items()]
            int_data=apply_colors(int_data)
            save_json(dyn_id,traj_id,int_type,int_data,resultspath)
##################


options = parser.parse_args()
overwrite=options.overwrite
dyn_dict_path=options.dyn_dict_path
dynamics_id=options.dynamics_id
ignore_publication=options.ignore_publication

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
gcdatapath_base="/protwis/sites/files/Precomputed/covid19/get_contacts_files/dynfreq_results/"
dynpath="/protwis/sites/files/Covid19Dynamics/"
resultspath="/protwis/sites/files/Precomputed/covid19/flare_plot/"

tot=len(dyn_dict.keys())
i=0
for dyn in sorted(dyn_dict.values(),key=lambda x:x["dyn_id"]):
    dyn_id=dyn["dyn_id"]
    #gcdata_path=os.path.join(gcdatapath_base,"dyn%(dynid)s/dyn%(dynid)s_dynamic.tsv" % {"dynid" :dyn_id} )
    gcdata_path_dyn=os.path.join(gcdatapath_base,"dyn%(dynid)s/" % {"dynid" :dyn_id} )
    if not os.path.isdir(gcdata_path_dyn):
        i+=1
        print("Skipping - dyn %s - Data not found" % dyn_id)
        continue
    print("dyn %s - %.1f%%"%(dyn_id , (i/tot)*100) )
    traj_files=dyn["files"]["traj"]
    dynfiles_traj=[]
    for mytrajd in traj_files:
        traj_id=mytrajd["id"]
        ex="%s_trj_%s_hb.json" % (traj_id,dyn_id)
        resultspath_hb=os.path.join(resultspath,"hb")
        if not overwrite and os.path.isfile(os.path.join(resultspath_hb,ex)):
            print("Skipping - dyn %s traj %s - already parsed"%(dyn_id,traj_id))
            continue
        traj_path=mytrajd["path"]
        if not os.path.isfile(traj_path):
            print("Skipping - dyn %s traj %s - Traj. file not found"%(dyn_id,traj_id))
            continue
        gcdata_path_dyn_traj=os.path.join(gcdata_path_dyn,"dyn%(dynid)s-%(trajid)s_dynamic.tsv" % {"dynid" :dyn_id, "trajid":traj_id} )
        if not os.path.isfile(gcdata_path_dyn_traj):
            print("Skipping - dyn %s traj %s - GetContacts results file not found"%(dyn_id,traj_id))
            continue
        #t=md.open(traj_path)
        #n_frames=t.__len__()
        #del t
        #gc.collect()
        #mytrajd["n_frames"]=n_frames
        dynfiles_traj.append(mytrajd)
    if dynfiles_traj:
        try:
            create_p_jsons(dynfiles_traj,gcdata_path_dyn,dyn_id,resultspath)
        except Exception as e:
            print("Skipping due to error: ",e)
    else:
        print("Skipping - Missing traj info")
    i+=1


