#!/usr/bin/python3
import sys,os,shutil,subprocess
import numpy
import pandas as pd

sys.path.insert(1,os.path.dirname(shutil.which('xtract')))
import edirect
#write a function for esearch and efetch
def search () :
  PN=input("Please enter the protein name:") or "glucose-6-phosphatase"
  OR=input("Please enter the organism name:") or "Aves"
  file=open("data.fa",'w')
  command="esearch -db protein -query "+"\""+PN+"[title] AND "+OR+"[Organism]"+"\""+" | efetch -format fasta"
  print(command)
  file.write(edirect.pipeline(command))
  file.close()
  sz=os.path.getsize("./data.fa")
  if sz == 0:
    return False
  else :
    return True

#check whether the data is null file
if search() == False : # the case that file is null
  choice=input("The output data is null.Do you want to change the key words and search again?(yes/no)") or "no" #ask user search again or not
  if choice == "yes":
     search()
  else :
     print("Quit. Go to help manual for help.")
     sys.exit()
else:
#the case the file is not null
  file=open("data.fa")
  data=file.read()
  seq=data.split(">")[1:] #split into seqs by ">"
  file.close()
  counts_seq=data.count(">") #calculate sequnces
  print("There're "+str(counts_seq)+" protein sequences have been found.")
  if counts_seq > 1000 :  # ask user whether use above 1000 seqs or not
     choice=input("There're more than 1000 sequences in the dataset, do you want to use only first 1000 sequences in the dataset?(yes/no)") or "yes" 
     if choice == "yes" :
        seq=seq[0:1000] #use first 1000 seqs
#find out how many species are included
  species_l=[] #create a list to save species
  for i in seq :
     if i.find('[') == -1 : #in case there's no species information
        continue 
     else:
        species=i.split("[")[1].split("]")[0] #grep species name
        species_l.append(species) #add it to the list
  species_num=len(set(species_l)) #count species without duplicates
  print("There're "+str(species_num)+" different specises in the dataset.")
  choice=input("Show the species list?(yes/no)") or "yes"
  if choice == "yes" :
     print("\n".join(set(species_l))) #show species list
  
#prosite
export_command="export EMBOSS_DATA=/localdisk/home/software/EMBOSS-6.6.0/share/EMBOSS/data/"
patmatmotifs_command=f"{export_command} && patmatmotifs -sprotein1 tem_input.txt -rformat2 nametable -outfile tem_motif.txt"
for i in seq :
  with open("tem_input.txt","w") as file :
    file.write(">"+i)
  process=subprocess.Popen(patmatmotifs_command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  stdout,stderr=process.communicate()
  print(stdout.decode())  # 打印标准输出结果
  print(stderr.decode())  # 打印错误输出结果
  with open("motif.out",'a') as file1 :
    file_size = os.path.getsize("motif.out")
    if file_size == 0 :
       with open("tem_motif.txt") as file2 :
         for line in file2 :
           if line.strip() and not line.startswith("#") :
              motif_line=[]
              motif_line=line.strip('\n').split()
              motif_line.append('\n')
              file1.write("\t".join(motif_line))
    with open("tem_motif.txt") as file2 :
       for line in file2 :
         if line.strip() and not line.startswith("#") and not line.startswith("USA") :
            motif_line=[]
            motif_line=line.strip('\n').split()
            motif_line.append('\n')
            file1.write("\t".join(motif_line))

#pandas show motifs associated
motif_df=pd.read_csv("motif.out",sep='\t',na_values=['-'])
Motif=motif_df['Motif'].value_counts()
with open("motifs_associated.out",'w') as file:
  file.write(Motif.to_string())
if len(Motif)==0 :
  print("Didn't find any motifs associated with these sequences from PROSITE.")
else :
  print(str(len(Motif))+" motifs are associated to sequence set.")
  choice=input("Show them now?(yes/no)") or "yes"
  if choice == "yes" :
     print(Motif.to_string())
