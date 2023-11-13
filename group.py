#!/usr/bin/python3
import sys,os,shutil,subprocess
import numpy
import pandas as pd

#put infomation on alignments
infoalign_command="infoalign -sequence Clusout.fa -only -name -idcount -change -outfile clusout.infoalign"
process = subprocess.Popen(infoalign_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()
#sort the Identity
column_names=["Name","Ident","%Change"]
df=pd.read_csv("clusout.infoalign",sep="\t",na_values=['-'],names=column_names)
df.sort_values('Ident',ascending=False,inplace=True)
max=str(df['Ident'].max())
min=str(df['Ident'].min())
median=str(df['Ident'].median())
print("The highest identity is "+max+" ,lowest is "+min +", and the median is "+median)

choice=input("Do you want to use a subset dataset filtered by Identities?(yes/no)") or "yes"

def Inputthreshold():
  set_h=input(f"Please enter the upper threshold value (<= {max}):") or max
  set_l=input(f"Please enter the lower threshold value (>= {min}):") or min
  return set_h,set_l

if choice == "yes" :
  check=True
  while(check):
    set_h, set_l = Inputthreshold()
    check=False
    try :
      set_h=float(set_h)
      set_l=float(set_l)
    except ValueError:
      print("the threshold is not a value. Please enter again.")
      check=True
  subset_df = df[df['Ident'].between(set_l, set_h)]
  l=subset_df.shape[0]
  print("There are "+str(l)+" sequences in the subset."+"\n"+"Now creating subset files...")
#creating subset
  file_data_origin=open("data.fa")
  data=file_data_origin.read()
  seq=data.split(">")[1:] #split into seqs by ">"
  file_data_origin.close()
  file_c=open("Clusout.fa")
  data2=file_c.read()
  seqc=data2.split(">")[1:]
  file_c.close()
  os.remove("subdata.fa")
  os.remove("subClusout.fa") 
  for index,name in subset_df['Name'].items() :
     with open("subdata.fa",'a') as file :
        subseq=seq[index]
        file.write(">"+subseq)
     with open("subClusout.fa",'a') as file :
        subseq2=seqc[index]
        file.write(">"+subseq2)
