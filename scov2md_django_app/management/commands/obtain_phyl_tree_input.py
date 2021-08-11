from django.core.management.base import BaseCommand, CommandError
import pickle
import os
import re
import seaborn as sns
import datetime
import pandas as pd

class Command(BaseCommand):
    help = "Transforms newick format to dictionary for the homepage phyl tree."
    def add_arguments(self, parser):
        parser.add_argument(
           '-i',
            dest='input_nwk',
            action='store',
            default='/protwis/sites/files/Covid19Data/Data/gisaid_timetree.nwk',
            help='Path to input file, which is the newick file'
        )
        parser.add_argument(
           '-m',
            dest='input_metadata',
            action='store',
            default='/protwis/sites/files/Covid19Data/Data/gisaid_metadata.tsv',
            help='Path to metadata file in TSV format'
        )
    def handle(self, *args, **options):

        def fix_metadata(strain_data):
            none_list=["","unknown","unknowne","unkown","?","not applicable"]
            sex_val=strain_data["Sex"]
            if sex_val.isnumeric():
                if strain_data["Age"] in none_list:
                    strain_data["Age"]=float(sex_val)
                strain_data["Sex"]=""
            elif sex_val=="Woman":
                strain_data["Sex"]="Female"
            for (datalabel,dataval) in strain_data.items():
                if type(dataval)==str and dataval.lower() in none_list:
                    strain_data[datalabel]=""
                if datalabel=="S1 mutations":
                    if dataval:
                        strain_data[datalabel]=str(dataval)
                    else:
                        strain_data[datalabel]="0"

            return strain_data

        def nwk_to_dict(newick,metadata):
            newick=newick.replace(" ","").replace("'","")
                                 #"(B:6.0,(A:5.0,C:3.0,E:-4.0)Ancestor1:5.0,D:11.0);"
            tokens = re.findall(r"([^:;,()\s]*)(?:\s*:\s*([-\d.e]+)\s*)?([,);])|(\S)", newick+";")
            def recurse(nextid = 0, parentid = -1): # one node
                thisid = nextid;
                children = []

                name, length, delim, ch = tokens.pop(0)
                genbank_accession=""
                date=""
                country=""
                if ch == "(":
                    while ch in "(,":
                        node, ch, nextid = recurse(nextid+1, thisid)
                        children.append(node)
                    name, length, delim, ch = tokens.pop(0)
                if length:
                    if float(length)<0:
                        length=0
                    if float(length)>1000:
                        length=0

                node_d={"id": thisid, "name": name, "length": float(length) if length else '', 
                        "parentid": parentid, "children": children , "is_terminal": False if children else True }

                if name:
                    if name in metadata.index: #This is for GSAID trees
                        strain_data=metadata.loc[name]
                        strain_data=fix_metadata(strain_data)
                        for (datalabel,dataval) in strain_data.items():
                            datalabel=datalabel.lower().replace(" ","_")
                            if datalabel=="collection_data":
                                datalabel="date"
                            node_d[datalabel]=dataval
                    elif  "||" in name:#this if for NCBI trees
                        (genbank_accession,x,country,state,date)=name.split("|")
                        node_d["genbank_accession"]=genbank_accession
                        node_d["country"]=country
                        node_d["date"]=date
                return node_d, delim, nextid

            return recurse()[0]       



        def split_date_obj(date):
            datespl=date.split("-")
            (year,month)=datespl[:2]
            if len(datespl)>2:
                day=datespl[2]
            else:
                day=1
            dateobj = datetime.datetime(int(year),int(month),int(day))
            return dateobj

        def get_all_values(el,val_list,all_values={}):
            for myvar in val_list:
                if myvar=="collection_data":
                    myvar="date"
                if myvar not in all_values:
                    all_values[myvar]=set()
                el_val=el.get(myvar)
                if el_val:
                    all_values[myvar].add(el_val)
            if el["children"]:
                for e in el["children"]:
                    all_values=get_all_values(e,val_list,all_values)
            return all_values


        def get_max_min_date(date_li):
            mindate=None
            maxdate=None
            for mydatestr in date_li:
                dateobj=split_date_obj(mydatestr)
                if (not mindate) or dateobj < mindate:
                    mindate=dateobj
                if (not maxdate) or dateobj > maxdate:
                    maxdate=dateobj
            return (mindate,maxdate)

        def date_to_num(date,mindate):
            if not date:
                return ""
            dateobj=split_date_obj(date)
            dateval=(dateobj - mindate).days
            if dateval <0:
                raise Exception("Error obtaining date value")
            return dateval

        def incorporate_datevalue_extnode(el,mindate,colormap_date):
            if el.get("date"):
                mydate=el["date"]
                mydate=mydate.replace("'","")
                dateval=date_to_num(mydate,mindate)
                el["datevalue"]=dateval
                el["color_date"]=colormap_date[dateval]
            if el["children"]:
                for e in el["children"]:
                    incorporate_datevalue_extnode(e,mindate,colormap_date)



 
        def get_color_schemes(phyl_dict):
            val_list=['age', 'gisaid_clade', 'clade', 'country', 'country_of_exposure', 'admin_division', 'division_of_exposure', 'host', 'originating_lab', 'pango_lineage', 'region', 'sex', 'submitting_lab', 'collection_data', 'author', 'region_of_exposure',"s1_mutations", "emerging_lineage", ]
            all_values=get_all_values(phyl_dict,val_list)
            colors_dict={}
            for (varname,varvals) in all_values.items():
                if varname =="date":
                    (mindate,maxdate)=get_max_min_date(varvals)
                    date_diff=(maxdate - mindate).days +1
                    colormap_date=sns.color_palette("Spectral",n_colors=date_diff).as_hex()
                    incorporate_datevalue_extnode(phyl_dict,mindate,colormap_date)
                    date_colors={}
                    date_colorsli=[]
                    for i,mydate in enumerate(sorted(varvals)):
                        mydate=mydate.replace("'","")
                        dateval=date_to_num(mydate,mindate)
                        date_colors[mydate]=colormap_date[dateval]
                        date_colorsli.append({"color":colormap_date[dateval],"value":mydate,"pos":i})
                    colors_dict[varname]=date_colors
                    colors_dict["date_list"]=date_colorsli
                    
                if varname=="age":
                    pre_age_groups=list(range(0,int(max(varvals)),10))
                    age_groups=[(age,age+9) for age in pre_age_groups]
                    colormap_age=sns.color_palette("inferno_r",n_colors=len(age_groups)).as_hex()
                    all_ages=list(range(0,int(max(varvals))+1))
                    age_colors={}
                    for myage in all_ages:
                        for i,(minage,maxage) in enumerate(age_groups):
                            if myage>= minage and myage<=maxage:
                                age_colors[myage]=colormap_age[i]
                    colors_dict[varname]=age_colors
                else:
                    any_not_number=[e for e in varvals if not str(e).isnumeric()]
                    if any_not_number:
                        varvals=sorted(varvals)
                    else:
                        varvals=sorted(varvals,key=lambda x:float(x))
                    mycolormap=sns.color_palette("Spectral",n_colors=len(varvals)).as_hex()
                    colors_dict[varname]={myval:mycolormap[i] for i,myval in enumerate(varvals)}
            return colors_dict


       #TO DO: automatize download of input

        input_nwk=options["input_nwk"]
        input_metadata=options["input_metadata"]
        if not os.path.isfile(input_nwk):
            raise Exception("Input file not found")
        if input_metadata:
            if not os.path.isfile(input_metadata):
                raise Exception("Metadata file not found")
            #metadata=pd.DataFrame.from_csv(input_metadata, sep='\t').dropna(how="all").fillna(value="")
            metadata=pd.read_csv(input_metadata, index_col=0, error_bad_lines=False,sep='\t',quoting=3).dropna(how="all").fillna(value="")
            #metadata=pd.read_table(input_metadata, index_col=0, error_bad_lines=False,sep='\t',quoting=3).dropna(how="all").fillna(value="")
        else:
            print("No metadata file provided.")
            metadata={}



        print("Generating plot...")
        with open(input_nwk, 'r') as content_file:
            nwk_str = content_file.read()
            phyl_dict=nwk_to_dict(nwk_str,metadata)

        #Obtain color scheme of all variables in the metadata
        colors_dict=get_color_schemes(phyl_dict)

        out_path="/protwis/sites/files/Covid19Data/Data/tree.data"
        print("Saving %s"% out_path)
        with open(out_path,"wb") as out_fileh:
            pickle.dump(phyl_dict,out_fileh)
        self.stdout.write(self.style.NOTICE("File successfully generated."))

        out_path_colors="/protwis/sites/files/Covid19Data/Data/colorscales.data"
        print("Saving %s"% out_path_colors)
        with open(out_path_colors,"wb") as out_fileh:
            pickle.dump(colors_dict,out_fileh)