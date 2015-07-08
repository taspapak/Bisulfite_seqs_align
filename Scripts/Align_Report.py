#!/usr/bin/python
from __future__ import division
from prettytable import PrettyTable
import operator
import re
import os

#Author:Papakonstantinou Anastasios

#Read "Total_Organisms" file(that is parsed) and modify in order to get all the organisms
def parse_organ():

    global spec

    l1=open("Scripts/Total_Organisms")
    lines1=l1.readlines()

    total_species=[]   

    for line in lines1[4:]:
       total_species.append(line+"\n")

    total_species = [w.replace('\n', '') for w in total_species]

    #Replace all Human occurences for Homo sapiens : Regex matching
    total_species=[re.sub(r'Human.*', "Homo sapiens", i) for i in total_species]
    total_species=[re.sub(r'H.sapiens.*', "Homo sapiens", i) for i in total_species]
    total_species=[re.sub(r'Homo.*', "Homo sapiens", i) for i in total_species]
    total_species=[re.sub(r'H.Sapiens.*', "Homo sapiens", i) for i in total_species]
    
    
    #Replace all Mouse occurences for Mus Musculus : Regex matching
    total_species=[re.sub(r'Mouse.*', "Mus musculus", i) for i in total_species]
    total_species=[re.sub(r'M.musculus.*', "Mus musculus", i) for i in total_species] 
    total_species=[re.sub(r'Mus.*', "Mus musculus", i) for i in total_species]


    #Replace all Canis Familiaris  occurences for Canis Familiaris : Regex matching
    total_species=[re.sub(r'Canis familiaris.*', "Canis familiaris", i) for i in total_species]
    total_species=[re.sub(r'Canis Familiaris.*', "Canis familiaris", i) for i in total_species]


    my_set1=set(total_species)
     
    #Keep the lenght of the species for the parsing of table after
    spec = len(my_set1)

    return total_species,my_set1,lines1


#Read "FW_Organisms" file(that is parsed) in order to get the organisms from FW_Conv_nt
def fw_parse():

    l2=open("Scripts/FW_Organisms")
    lines2=l2.readlines()

    fw_species=[]

    for line in lines2[1:]:
       fw_species.append(line+"\n")

    fw_species = [w.replace('\n', '') for w in fw_species]

    #Replace all Human occurences for Homo sapiens : Regex matching
    fw_species=[re.sub(r'Human.*', "Homo sapiens", i) for i in fw_species]
    fw_species=[re.sub(r'H.sapiens.*', "Homo sapiens", i) for i in fw_species]
    fw_species=[re.sub(r'H.Sapiens.*', "Homo sapiens", i) for i in fw_species]
    fw_species=[re.sub(r'Homo.*', "Homo sapiens", i) for i in fw_species]

    #Replace all Mouse occurences for Mus Musculus : Regex matching
    fw_species=[re.sub(r'Mouse.*', "Mus musculus", i) for i in fw_species]
    fw_species=[re.sub(r'M.musculus.*', "Mus musculus", i) for i in fw_species]
    fw_species=[re.sub(r'Mus.*', "Mus musculus", i) for i in fw_species]


    #Replace all Canis Familiaris  occurences for Canis Familiaris : Regex matching
    fw_species=[re.sub(r'Canis familiaris.*', "Canis familiaris", i) for i in fw_species]
    fw_species=[re.sub(r'Canis Familiaris.*', "Canis familiaris", i) for i in fw_species]


    my_set2=set(fw_species)

    return fw_species,my_set2,lines2

#Read "RC_Organisms" file(that is parsed) in order to get the organisms from RC_Conv_nt
def rc_parse():

    l3=open("Scripts/RC_Organisms")
    lines3=l3.readlines()

    rc_species=[]

    for line in lines3[1:]:
       rc_species.append(line+"\n")

    rc_species = [w.replace('\n', '') for w in rc_species]

    #Replace all Human occurences for Homo sapiens : Regex matching
    rc_species=[re.sub(r'Human.*', "Homo sapiens", i) for i in rc_species]
    rc_species=[re.sub(r'H.sapiens.*', "Homo sapiens", i) for i in rc_species]
    rc_species=[re.sub(r'Homo.*', "Homo sapiens", i) for i in rc_species]
    rc_species=[re.sub(r'H.Sapiens.*', "Homo sapiens", i) for i in rc_species]
    

    #Replace all Mouse occurences for Mus Musculus : Regex matching
    rc_species=[re.sub(r'Mouse.*', "Mus musculus", i) for i in rc_species]
    rc_species=[re.sub(r'M.musculus.*', "Mus musculus", i) for i in rc_species]
    rc_species=[re.sub(r'Mus.*', "Mus musculus", i) for i in rc_species]
    
    #Replace all Canis Familiaris  occurences for Canis Familiaris : Regex matching
    rc_species=[re.sub(r'Canis familiaris.*', "Canis familiaris", i) for i in rc_species]
    rc_species=[re.sub(r'Canis Familiaris.*', "Canis familiaris", i) for i in rc_species]

    my_set3=set(rc_species)

    return rc_species,my_set3,lines3

#Prettytable organization
def table_elem(species,my_set,fw_species,my_set2,rc_species,my_set3):
    
    global spec

    l1=open("Scripts/Total_Organisms")
    lines1=l1.readlines()

    #Read the number of queries
    reads=lines1[0].split("R",1)[0].split(":",1)[1]

    tab = PrettyTable(["Organisms", "Hits", "Percentage(%)"])
    tab.padding_width = 1


    hits,fw_hits,rc_hits=0,0,0
    total=0
    
    my_dict={}
    fw_dict={}
    rc_dict={}

    fw_conv=[]
    rc_conv=[]

    perc=[]


    #############
    for i in my_set:
       indices = [j for j, x in enumerate(species) if x == i]
       my_dict.update({i: len(indices)}) 
       hits+= len(indices)

    for i in my_set2:
       fw_indices = [j for j, x in enumerate(fw_species) if x == i]
       fw_dict.update({i: len(fw_indices)}) 
       fw_hits+= len(fw_indices)

    for i in my_set3:
       rc_indices = [j for j, x in enumerate(rc_species) if x == i]
       rc_dict.update({i: len(rc_indices)}) 
       rc_hits+= len(rc_indices)
    

    ###########################
    for w in sorted(my_dict):

        perc.append(my_dict[w])
        #test2.append(round(my_dict[w]/hits*100,2))
        
        if w not in fw_dict:
           fw_dict.update({w:0})

        if w not in rc_dict:
           rc_dict.update({w:0})
           
        tab.add_row([w ,my_dict[w], round(my_dict[w]/int(reads)*100,2)])
        
        total += my_dict[w]/hits

    for s in sorted(fw_dict):
        fw_conv.append(fw_dict[s])
    
    for p in sorted(rc_dict):
        rc_conv.append(rc_dict[p])
    

    unconv=[sum(x) for x in zip(fw_conv, rc_conv)]

    un_conv=[i - j for i, j in zip(perc, unconv)]

    un_hits=sum(un_conv)

    fw_perc=[round(i / j*100,2) for i, j in zip(fw_conv,perc)]
    rc_perc=[round(i / j*100,2) for i, j in zip(rc_conv,perc)]
    un_perc=[round(i / j*100,2) for i, j in zip(un_conv,perc)]

    fw_conv2=map(str,fw_conv)
    fw_perc2=map(str,fw_perc)


    fw_perc3=[]


    for n in fw_perc2:
        fw_perc3.append("(" + n + "%)")

 
    c = map(" ".join, zip(fw_conv2, fw_perc3))
    

    rc_conv2=map(str,rc_conv)
    rc_perc2=map(str,rc_perc)

    rc_perc3=[]


    for n in rc_perc2:
        rc_perc3.append("(" + n + "%)")


    c2 = map(" ".join, zip(rc_conv2, rc_perc3))

    un_conv2=map(str,un_conv)
    un_perc2=map(str,un_perc)

    un_perc3=[]


    for n in un_perc2:
        un_perc3.append("(" + n + "%)")


    c3 = map(" ".join, zip(un_conv2, un_perc3))


    tab.add_column("FW_Conv_nt",c)
    tab.add_column("RC_Conv_nt",c2)
    tab.add_column("Unconv_nt" ,c3)

    ##########################################################

    #Alternative output:Percentage of each database as factor of the percentage of the organism

    #fw_perc2=[round(i * j/100,2) for i, j in zip(fw_perc,test2)]
    #rc_perc2=[round(i * j/100,2) for i, j in zip(rc_perc,test2)]
    #un_perc2=[round(i * j/100,2) for i, j in zip(un_perc,test2)]
    
    no_hits=int(reads)-hits

    
    tab.add_row(["#","#","#","#","#","#"])
    tab.add_row(["Total Queries:" ,reads,round(hits/int(hits)*100,2),"-","-","-"])
    tab.add_row(["No hits",str(no_hits),round(no_hits/int(reads)*100,2),"-","-","-"])
    tab.add_row(["Hits",str(hits), round(hits/int(reads)*100,2) ,str(round(fw_hits/hits*100,2))+ "%",str(round(rc_hits/hits*100,2)) + "%",str(round(un_hits/hits*100,2)) + "%"])
    
    
    long_org=max(my_set,key=len)
    long_hit=str(max(my_dict.iteritems(),key=operator.itemgetter(1))[1])
    long_fw=max(c,key=len)
    long_rc=max(c2,key=len)
    long_un=max(c3,key=len)
    
       
    bounds1=8*(len(long_org)+4+len(long_hit)+4+16+len(long_fw)+5+len(long_rc)+3+len(long_un)+3)
    bounds2=3*(len(long_org)+4+len(long_hit)+4+16+len(long_fw)+3+len(long_rc)+3+len(long_un)+3) 

    fin_tab=tab.get_string(start=0,end=spec,sortby="Hits",reversesort=True)
    fin_tab2=tab.get_string(start=spec,end=spec+4,sortby="Organisms",reversesort=True)


    return fin_tab2[0:bounds1],fin_tab[bounds2:len(fin_tab)]


#Write Report file
def write_rep(lines,table,table2):
    
    with open('Report.txt', 'w') as f2:
        
         f2.write("\n                                       !Organisms Alignment Report!\n\n")
    
         for line in lines[0:2]:
            f2.write("#" + line)
         
         f2.write(str(table)+"\n")
         f2.write(str(table2)+"\n")

    f2.close()
    
    #os.system("rm Total_Organisms FW_Organisms RC_Organisms")
    

#Methods call
def main():

    organisms,unique,lines=parse_organ()
    
    fw_organisms, fw_unique , fw_lines = fw_parse()
    rc_organisms, rc_unique , rc_lines = rc_parse()

    prettytab,prettytab2=table_elem(organisms,unique,fw_organisms,fw_unique,rc_organisms,rc_unique)
    
    write_rep(lines,prettytab,prettytab2)


if __name__ == "__main__":
    main()

