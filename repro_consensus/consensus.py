
import csv
import glob
import numpy as np



#Reading the OG/gold standard results
with open('./repro_consensus/results/lineage_report.csv', newline='') as csvfile:
    true_lineage_report = next(csv.DictReader(csvfile))

with open('./repro_consensus/results/annotations.tsv', newline='') as csvfile:
    true_annotations = next(csv.DictReader(csvfile, delimiter="\t"))


#Reading the groups
groups = glob.glob(f'./repro_consensus/results_participants/*', recursive=False)

print(f"There are {len(groups)} groups who have sent their results to participate in the reproducibility consensus.\n")



print("Annotation Consensus")
print("--------------------")
indentical, indentical_clade = [], []
for g in groups:
    csv_file = glob.glob(f'{g}/*.tsv', recursive=False)
    with open(csv_file[0], newline='') as csvfile:
        g_annotations = next(csv.DictReader(csvfile, delimiter="\t"))
   
    #Check if they are totally identical
    if(g_annotations==true_annotations):
        indentical.append(1)
    else:
        indentical.append(0)
   
    #Identical clade
    if(g_annotations["clade"]==true_annotations["clade"]):
        indentical_clade.append(1)
    else:
        indentical_clade.append(0)
    
print(f"There are {len(indentical)} groups ({len(indentical)/len(groups)*100:.0f}%) which have sent in the annotations file")
print(f"There are {np.sum(indentical)} groups ({np.sum(indentical)/len(indentical)*100:.0f}%) which have exactly the same annotation file to the gold standard")
print(f"There are {np.sum(indentical_clade)} groups ({np.sum(indentical_clade)/len(indentical)*100:.0f}%) which have exactly the same biological interpretation (clade) to the gold standard")
print()

print("Lineage Report Consensus")
print("--------------------")
indentical, indentical_lineage = [], []
for g in groups:
    csv_file = glob.glob(f'{g}/*.csv', recursive=False)
    with open(csv_file[0], newline='') as csvfile:
        g_lineage_report = next(csv.DictReader(csvfile))
   
    #Check if they are totally identical
    if(g_lineage_report==true_lineage_report):
        indentical.append(1)
    else:
        indentical.append(0)
   
    #Identical lineage
    if(g_lineage_report["lineage"]==true_lineage_report["lineage"]):
        indentical_lineage.append(1)
    else:
        indentical_lineage.append(0)
    
print(f"There are {len(indentical)} groups ({len(indentical)/len(groups)*100:.0f}%) which have sent in the lineage report")
print(f"There are {np.sum(indentical)} groups ({np.sum(indentical)/len(indentical)*100:.0f}%) which have exactly the same lineage report to the gold standard")
print(f"There are {np.sum(indentical_lineage)} groups ({np.sum(indentical_lineage)/len(indentical)*100:.0f}%) which have exactly the same biological interpretation (lineage) to the gold standard")




#print(true_lineage_report)
#print(true_annotations)
