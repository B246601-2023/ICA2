#!/usr/bin/python3
import os,sys,subprocess,shutil,random,re
import pandas as pd
sys.path.insert(1,os.path.dirname(shutil.which('xtract')))
import edirect



answers={}#create a dictionary to contain inputs
#write a function for esearch and efetch
def search () : #function
    answers["PN"]=input("Please enter the protein name : ") or "glucose-6-phosphatase"
    answers["field1"]=input(f"Please eneter the search field for \"{answers['PN']}\" [default: title] : ") or "title"
    answers["OR"]=input("Please enter the organism name : ") or "Aves"
    file=open("data.fa",'w')
    command="esearch -db protein -query "+"\""+answers["PN"]+f"[{answers['field1']}] AND "+answers["OR"]+"[Organism]"+"\""+" | efetch -format fasta"
    print(command) # show the command to usrs
    file.write(edirect.pipeline(command)) # put sequences got to "data.fa"
    file.close()
    sz=os.path.getsize("./data.fa") # determine whether got any sequences
    if sz == 0:
        return False # to give an indicator to start question loop
    else :
        return True
def create_folder(fpath): # create folder if not exists
    if not os.path.exists(fpath):
        os.makedirs(fpath)
        print(f"Folder '{fpath}' created.\n")
    return 0

def clear_file(fpath): # clear file if exists
    if os.path.exists(fpath):
        os.remove(fpath)
    return 0

def clear_folder(folder_path): # clear contents in a folder
    if os.path.exists(folder_path): # determine whether the directory exists
        for file_name in os.listdir(folder_path): # iterate through all files and subfolders within a folder
            file_path = os.path.join(folder_path, file_name)
            if os.path.isfile(file_path): #determin file or directory
                os.remove(file_path) #If it is a file, delete the file
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path) #If it is a directory, delete the directory and its contents recursively
    return 0
def Inputthreshold(max,min): # get thresholds from usr
  set_h=input(f"Please enter the upper threshold value (<= {max}) : ") or max
  set_l=input(f"Please enter the lower threshold value (>= {min}) : ") or min
  return set_h,set_l

def create_infoalign(): # create infoalign file for sequences
    infoalign_command = "infoalign -sequence Clusout.fa -only -name -idcount -change -outfile clusout.infoalign" # only need name idcount and change columns
    infoalign = subprocess.Popen(infoalign_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = infoalign.communicate()
    print(stdout.decode())
    print(stderr.decode())

    column_names = ["Name", "Ident", "%Change"] # create colnames for the table
    df = pd.read_csv("clusout.infoalign", sep="\s+", na_values=['-'], names=column_names)  # load dataframe
    #the clusout.infoalign are not in the perfect format use '\t' to seperate, so use "\s+" here.
    return df # return a dataframe

def creating_seq(path): # load sequences into a list without ">"; an element stands for a seq
    file_data = open(path)
    data = file_data.read()
    seq = data.split(">")[1:]  # split into seqs by "> ; skip the first '' element
    file_data.close()
    return seq # return a list

def search_add_seq_onebyone(name,seqset,path): # a fucntion to check if a accession id can be found ; use in choice-type2 ;
    # name : accession id; seqset : a sequence list to be checked ; path : distinguish original sequecnes or aligned seqeunces
    found = False # an indicator
    for seq in seqset: #iterate all sequences
        if re.search(name, seq): # if id exists, write relating sequences into a susbet file
            print("One sequence found. Adding to the subset now ...\n")
            found = True
            with open(f"./subsets/{path}{i}.fa", 'a') as file: # i indicates current subset number
                subseq = seq
                file.write(">" + subseq)
    if found == False: # if id not found
        print("No sequences found.\n")
    return 0


def replace_ambiguous (filepath): # replace X U B in sequences
    with open(filepath) as file:
        lines = file.readlines()
    processed_sequences = [] # List to store processed sequences
    for line in lines: # process each line
        if line.startswith('>'):
            # If the line starts with '>', it indicates the beginning of a new sequence
            processed_sequences.append(line) # Save the previous sequence (if any) and start processing a new sequence
        else:
            # For lines that do not start with '>', replace ambiguous characters 'X', 'U', 'B' with '-'
            current_sequence = line.replace('X', '-').replace('U', '-').replace('B', '-')
            processed_sequences.append(current_sequence)

    # save results to a new file
    with open("tem_no.fasta", 'w') as output_file:
        output_file.writelines(processed_sequences)


def choi_type(i): # a function to create subsets with 3 method ; i is the number of current subset being processed
    print(f"Now creating subset {i}...\n")
    choi_type = input("How to create the subset?\n\
    1.Filter by sequence identity\n\
    2.Enter sequence ACCESSION by hand\n\
    3.Choose n sequences by random\n\
    Please enter 1/2/3 to choose : ")
    if choi_type == "1":
        if not os.path.exists("clusout.infoalign") :
            df = create_infoalign() # create dataframe with identity scores
        else:
            column_names = ["Name", "Ident", "%Change"]
            df = pd.read_csv("clusout.infoalign", sep="\s+", na_values=['-'], names=column_names)
        df.sort_values('Ident', ascending=False, inplace=True)  # sort by Ident
        max = str(df['Ident'].max())  # get the max value of Ident
        min = str(df['Ident'].min())  # get the min value of Ident
        median = str(df['Ident'].median())
        print(
            "The highest identity is " + max + " ,lowest is " + min + ", and the median is " + median+"\n"
        )  # tell usr the range of Ident

        # ensure usr put in correct thresholds
        check = True
        while (check):
            set_h, set_l = Inputthreshold(max,min)
            check = False
            try:
                set_h = float(set_h)
                set_l = float(set_l)
            except ValueError: # if value error keep loop
                print("the threshold is not a value. Please enter again.\n")
                check = True
        subset_df = df[df['Ident'].between(set_l, set_h)] # select sequences meet requirments
        l = subset_df.shape[0] # count how many sequences are chosen
        if l == 0 :
            print("No sequences found in this subset...\n")
            nextstep = input("Continue(0) or create again(2)? : ") or "0" # if no sequence, recreate or continue.
            return nextstep,i
        print("There are " + str(l) + " sequences in the subset.\n" + "Now creating subset files...\n")
        create_folder("./subsets")
        seq_d = creating_seq("data.fa") # get seq list from data.fa
        seq_c = creating_seq("Clusout.fa")# get seq list from Clusout.fa
        clear_file(f"./subsets/subdata{i}.fa")
        clear_file(f"./subsets/subClusout{i}.fa")
        for index, name in subset_df['Name'].items(): # index in df is the same in seq list
            with open(f"./subsets/subdata{i}.fa", 'a') as file:
                subseq = seq_d[index]  # get the correct seq and write into a susbset file
                file.write(">" + subseq)
            with open(f"./subsets/subClusout{i}.fa", 'a') as file:
                subseq2 = seq_c[index]
                file.write(">" + subseq2)
        nextstep=input(f"Subset {i} successfully created! Continue to next function(0) or create another subset(1)? : ") or "0"
        i=i+1
    elif choi_type == "2":
        if not os.path.exists("clusout.infoalign") : # create clustout.infoalign
            df = create_infoalign()
        else:
            column_names = ["Name", "Ident", "%Change"]
            df = pd.read_csv("clusout.infoalign", sep="\s+", na_values=['-'], names=column_names)
        show_choice = input("Show all sequences ACCESSION first ? (yes/no) : ") or "no"
        if show_choice == "yes" :
            print(df['Name'].to_string(index=False))
        create_folder("./subsets")
        seq_d = creating_seq("data.fa")
        seq_c = creating_seq("Clusout.fa")
        clear_file(f"./subsets/subdata{i}.fa")
        clear_file(f"./subsets/subClusout{i}.fa")
        while True:  # keep asking until usr input 'q'
            name = input("Please enter the ACCESSION (ONE at a time!; 'q' to quit) : ") or "q"
            if name == "q" :
                print("Input complete.\n") # leave the loop
                break
            search_add_seq_onebyone(name,seq_d,path="subdata")
            search_add_seq_onebyone(name,seq_c,path="subClusout")
        if not os.path.exists(f"./subsets/subdata{i}.fa") or not os.path.exists(f"./subsets/subClusout{i}.fa") : # if no sequences are put into the current dataset
            print("No sequence in current subset.")
            nextstep = input("Continue to next function (0) or create this subset again(2) : ?") or "0"
            return nextstep,i # i stay with the same subset
        nextstep = input(f"Subset {i} successfully created! Continue(0) or create another subset(1) : ?") or "0"
        i=i+1 # move to next subset
    elif choi_type == "3":
        if not os.path.exists("clusout.infoalign") :
            df = create_infoalign()
        else:
            column_names = ["Name", "Ident", "%Change"]
            df = pd.read_csv("clusout.infoalign", sep="\s+", na_values=['-'], names=column_names)
        l = df.shape[0] # get number of sequences
       # ensure the usr give the correct sequence amount
        check = True
        while (check): # stay in loop unless n within (0,number of sequences)
            n = input(f"Please enter the number of sequences to grab (n <= {l}) : ") or "10"
            check = False
            try:
                if int(n) > l or int(n) <= 0:
                    n = "null"
                else :
                    n = int(n)
                n = int(n)
            except ValueError:
                print("n is not valid. Please enter again.\n")
                check = True
        uniq_num = random.sample(range(0, l), n) # generate a random number set
        create_folder("./subsets")
        seq_d = creating_seq("data.fa") # get seq list from data.fa
        seq_c = creating_seq("Clusout.fa") # get seq list from Clusout.fa
        clear_file(f"./subsets/subdata{i}.fa")
        clear_file(f"./subsets/subClusout{i}.fa")
        for index in uniq_num: # use numberset as index to get seqeunce ; iterate all numbers
            with open(f"./subsets/subdata{i}.fa", 'a') as file:
                subseq = seq_d[index]
                file.write(">" + subseq)
            with open(f"./subsets/subClusout{i}.fa", 'a') as file:
                subseq2 = seq_c[index]
                file.write(">" + subseq2)
        nextstep = input(f"Subset {i} successfully created! Continue to next function(0) or create another subset(1)? : ") or "0"
        i=i+1 # move to next subset
    else:
        nextstep=input("No subsets created. Continue to next function(0) or create another subset(1)? : ") or "0" # i stay with the same subset
    return nextstep,i

#plotcon
def plotco(winsize,path,p) : # p: give current plot a number
  if p == 1:
     title="Full"
  else:
     title=path
  plotcon_commands=[
                      f"plotcon {path} -winsize {winsize} -graph x11 -gsubtitle {title}", #I cannot prompt and save them at the same time.
                      f"plotcon {path} -winsize {winsize} -graph png -gdirectory ./out_plots/conservation_plot -goutfile plot{p} "
                   ]
  for plotcon_command in plotcon_commands :
     plotcon=subprocess.Popen(plotcon_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     stdout, stderr = plotcon.communicate()
     print(stdout.decode())
     print(stderr.decode())
  return 0

#prosite
def getmotifs(path,c) : # path : used to read sequences from ; c : decide which subset now is processing
    export_command="echo $EMBOSS_DATA"
    #export_command = "export EMBOSS_DATA=/localdisk/home/software/EMBOSS-6.6.0/share/EMBOSS/data/"
    patmatmotifs_command=f"{export_command} && patmatmotifs -sprotein1 tem_input.txt -rformat2 nametable -outfile tem_motif.txt" # get a motif output for one sequence
    seq = creating_seq(path) #create a list with all sequences
    clear_file("motif.out")
    for i in seq : # Iterate through all sequences, process them and save their results into one summary file with table format
        with open("tem_input.txt","w") as file :
            file.write(">"+i) # put current sequence into a file to process
        patmatmotif=subprocess.Popen(patmatmotifs_command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=patmatmotif.communicate() #execute
        print(stdout.decode())
        print(stderr.decode())
        with open("motif.out",'a') as file1 : # the summary file
            file_size = os.path.getsize("motif.out")
            if file_size == 0 : # create with table title first("USA Motif ...")
                with open("tem_motif.txt") as file2 : # open the motif result for one sequence
                    for line in file2 : # iterate all lines in tem_motif.txt
                        if line.strip() and not line.startswith("#") : # find lines with both title and contents
                            motif_line=[] #create a list
                            motif_line=line.strip('\n').split() # sperate lines in tem_motif.txt by'', and use a list to hold elements .
                            motif_line.append('\n') # add a '\n' between elements from each line
                            file1.write("\t".join(motif_line)) #join the elements with '\t' , form a table saved in file 'motif.out'
            with open("tem_motif.txt") as file2 : # after the title has been added to the table, only grab useful contents
                for line in file2 :
                    if line.strip() and not line.startswith("#") and not line.startswith("USA") :
                        motif_line=[]
                        motif_line=line.strip('\n').split()
                        motif_line.append('\n')
                        file1.write("\t".join(motif_line))

#pandas show motifs associated
    motif_df=pd.read_csv("motif.out",sep='\t',na_values=['-']) # read the summary table
    Motif=motif_df['Motif'].value_counts() # extract motif column
    create_folder("./outputfiles")
    with open(f"./outputfiles/motifs_associated{c}.out",'w') as file:
        file.write(Motif.to_string()) # write df into a file
    if len(Motif)==0 :
        print(f"Didn't find any motifs associated with dataset{c} from PROSITE.\n")
    else :
        print(str(len(Motif))+f" motifs are associated to dataset{c}.\n")
        choice=input("Show them now?(yes/no) : ") or "yes"
        if choice == "yes" :
             print(Motif.to_string()+"\n")
    clear_file("tem_input.txt") # clear temporary files
    clear_file("tem_motif.txt")



def sigcleave(filepath,output_filepath) :
    while True:
        weight=input("Minimum scoring weight value for the predicted cleavage site (from 0.00 to 100.000) [default:3.5] : ") or 3.5
        try: # ensure inputs are valid
            if float(weight) <= 100 or float(weight) >= 0 :
                break
            print("Please enter invalid value (from 0.00 to 100.000).\n")
        except ValueError:
            print("Please enter invalid value (from 0.00 to 100.000).\n")
    replace_ambiguous(filepath)
    sig_command=f"sigcleave tem_no.fasta -minweight {weight} -rformat nametable {output_filepath}"
    process = subprocess.Popen(sig_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout.decode())
    print(stderr.decode())
    return 0


#prettyplot:Draw a sequence alignment with pretty formatting
def prettyplot(path,subset): # path : the dataset to use; subset : number of the subset to use
    while True:
        ratio=input("Please specity a plurality ratio for a consensus match (default=0.5, 0 <= ratio <= 1) : ") or 0.5
        try :
            float(ratio)
            if float(ratio) <= 1 and float(ratio) >= 0 :
                break
            else :
                print("The ratio is invalid.\n")
        except ValueError:
            print("The ratio is invalid.\n")

    prettyplot_command=f"prettyplot {path} -ratio {ratio} -graph png -gdirectory ./out_plots/prettyplot_{subset}"
    process = subprocess.Popen(prettyplot_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    print(stdout.decode())  # print standard output
    print(stderr.decode())  # print error output
    return 0


#initialize
clear_folder("./outputfiles")
clear_file("clusout.infoalign")
clear_file("motif.out")
clear_file("Clusout.fa_10")
clear_folder("./subsets")
clear_folder("./out_plots")

#check whether the data is null file
nonull=search() # get the indicator
while nonull == False : # the case that file is null, keep loop until te file is not null
    choice=input("The output data is null.Do you want to change the key words and search again? (yes/no) : ") or "no" #ask user search again or not
    if choice == "yes":
        nonull=search()
    else :
        print("Quit the script. Go to help manual for help.")
        sys.exit() #exit script

# the case that file is not null
seq=creating_seq("data.fa")
counts_seq=len(seq) #calculate sequnces
print("There're "+str(counts_seq)+" protein sequences have been found.\n")
if counts_seq > 1000 :  # ask user whether use above 1000 seqs or not
    choice=input("There're more than 1000 sequences in the dataset, do you want to use only first 1000 sequences in the dataset? (yes/no) : ") or "yes"
    if choice == "yes" :
        seq=seq[0:1000] #use first 1000 seqs
        for s in seq : #change data.fa with ">"
            with open("data.fa",'w') as file :
                file.write(">"+s)
#find out how many species are included
species_l=[] #create a list to save species
for i in seq :
    if i.find('[') == -1 : # in case there's no species information then neglect
        continue
    else:
        species=i.split("[")[1].split("]")[0] # grep species name
        species_l.append(species) # add it to the list
species_num=len(set(species_l)) # count species without duplicates
print("There're "+str(species_num)+" different specises in the dataset.\n")
choice=input("Show the species list? (yes/no) : ") or "yes" #ask
if choice == "yes" :
    print("\n".join(set(species_l))+"\n") # show species list


#clustalo
clustalo_command = "clustalo -i data.fa --guidetree-out=guidetree.out -o Clusout.fa --outfmt=fasta -v --force"
print("Clustering...")
clustalo = subprocess.Popen(clustalo_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = clustalo.communicate() # execute clustalo command
print(stdout.decode())
print(stderr.decode())

#creating subset
while True : # ensure inputs valid
    answers["choi_subset"] = input("Do you want to create a subset from all sequences? (yes/no) : ") or "no"
    if answers["choi_subset"] == "no" :
        check_main = False # create an indicator for next loop
        break
    elif answers["choi_subset"] == "yes":
        check_main = True
        break
    else :
        print("Invalid input.\n")
i=1 # i is the number of the current subset
while (check_main):
    nextstep,i=choi_type(i) # create subset
    if nextstep == "0" :
        break # exit loop
    elif nextstep == "1" :
        continue # move to next loop
    elif nextstep == "2" :
        print(f"Recreating subset {i}...\n")
    else :
        print("Invalid answer. Please enter again.\n")


#plotcon
print("Starting plotcon function now...\n")
create_folder("./out_plots/conservation_plot") # create a dir for outputs
fullset="./Clusout.fa"

if answers["choi_subset"] == "yes" and os.path.exists("./subsets/subClusout1.fa"): # Subsets exists, ask usr to choose which dataset to use
    subset = [f for f in os.listdir("./subsets") if f.startswith("subClusout")]  # make a list of clustered subset files
    while True : # keep ask unless get the valid answer
        choice=input("Which dataset to use for generating plotcon? (full/subset/both) :" ) or "both"
        if choice in ["both","subset","full"] :
            break
        else:
            "Invalid choice. Enter again. \n"
else: # Subsets don't exist, use original dataset
    print("Use full set to generate plotcon.")
    choice="full" #use full dataset because there's no subset
winsize=input("Please enter the windowsize value [default=4] : ") or 4
if choice == "both" :
  plotco(winsize,fullset,1) # plotcon for full set
  p=2 # count+1
  for path in subset: # plotcon for each subset
      plotco(winsize,f"./subsets/{path}",p)
      p=p+1
elif choice == "subset" :
    for path in subset:
        p=2
        plotco(winsize, f"./subsets/{path}", p)
        p=p+1
elif choice == "full" :
  plotco(winsize,fullset,1)


#use getmotifs

print("Starting getmotifs function now...\n")
if answers["choi_subset"] == "no" or not os.path.exists("./subsets/subClusout1.fa"): # case subsets didn't exist
    path = "data.fa"
    getmotifs(path,c=1)
else :
    c=1 # start from subset 1
    file_list = sorted(os.listdir("./subsets")) # keep susbets on numeric sort
    for file in file_list :
        if re.search("subdata",file) : # put in original sequences not aligned one
            file_name=re.findall("subdata.+",file)[0]
            path = f"./subsets/{file_name}"
            getmotifs(path,c)
            c=c+1 # move to next subset
    print("Quit getmotifs function...\n")
    print(str(c)+" Motifs output has been put in ./outputfiles."+"\n")
clear_file("motif.out")

#use sigcleave
print("Starting sigcleave function...\n")
count=0 # the number of file processed
if answers["choi_subset"] == "yes" and os.path.exists("./subsets/subClusout1.fa"): #if the subset exists
    file_list = sorted(os.listdir("./subsets"))# keep susbets on numeric sort
    while True: # stay in loop until usr asked to quit
        print("Subset list : \n")
        print("\n".join(file_list)+"\n") #show dir list
        subset=input("Choose the subset to use to generate report on signal cleavage sites (Please enter the number 1/2/3...; 'q' to quit) : ")
        if subset=="q" :
            print("Quit sigcleave function...\n")
            break
        subset_path=f"./subsets/subClusout{subset}.fa"
        find_anyfile=False # indicator to show whether find the subset usr asked
        for file in file_list:
            if re.search(file,subset_path) :
                find_anyfile=True
                sigcleave(subset_path,f"./outputfiles/subset{subset}.sig") # generate report
                count=count+1 # move to next susbet
                clear_file("tem_no.fasta")
        if not find_anyfile :
            print("No such subset.\n")
else:
    choice = input("Skip this function? ('q' to quit ; Press anything to continue) : ") or "q"
    if choice != "q" :
        count = 1  # only one file processed
        subset_path = "Clusout.fa"
        sigcleave(subset_path, "./outputfiles/fullset.sig")
        clear_file("tem_no.fasta")
    else :
        count = 0

print(f"There're {count} reports saved in dir 'outputfiles'\n")
print("Quit sigcleave function...\n")

#use prettyplot
print("For better utilization of prettyplot, only subsets with less than 10 sequences will be used to generate graphs.\n")
choice=input("Would you like to regenerate the subset? (yes/no) : ") or "no"
if choice == "yes" : # the same as create subsets part
    while True:
        nextstep,i=choi_type(i)
        if nextstep == "0" :
            break
        elif nextstep == "1" :
            continue
        elif nextstep == "2" :
            print(f"Recreating subset {i}...\n")
        else :
            print("Invalid answer. Please enter again.\n")

print("Starting prettyplot function now...\n")
if answers["choi_subset"] == "yes" : #if the subset exists
    file_list = sorted(os.listdir("./subsets"))
    while True: # stay in loop unless usr aked to quit
        print("Subset list : \n")
        print("\n".join(file_list)+"\n")
        subset=input("Choose the subset to use to generate alignment plots (Please enter the number 1/2/3...; 'q' to quit) : ")
        if subset=="q" :
            print("Quit Prettyplot function...\n")
            break
        subset_path=f"./subsets/subClusout{subset}.fa"
        create_folder(f"./out_plots/prettyplot_{subset}")
        find_anyfile=False
        for file in file_list:
            if re.search(file,subset_path) : # search for subset usr asked
                seq=creating_seq(subset_path)
                find_anyfile=True
                if len(seq) >= 10 : # use first 10 seq if seq number larger than 10
                    print(f"More than 10 sequences in subest{subset}. Automatically use first 10 sequences...\n")
                    with open(subset_path+"_10",'w') as file:
                        seq=seq[0:10]
                        seq_write=[f">{item}" for item in seq] # add ">" to each sequence
                        file.write("".join(seq_write))
                    print(f"Generating alignment plots for subset{subset}...\n")
                    prettyplot(subset_path+"_10",subset)
                elif len(seq)<=10 and len(seq)>=0 :
                    print(f"Generating alignment plots for subset{subset}...\n")
                    prettyplot(subset_path,subset)
        if not find_anyfile :
            print("No such subset.\n")
else:
    choice=input("Skip this function? ('q' to quit ; Press anything to continue) : ") or "q"
    if choice != "q" :
        subset_path = "Clusout.fa"
        seq = creating_seq(subset_path)
        subset = "full"
        create_folder(f"./out_plots/prettyplot_{subset}")
        if len(seq) >= 10: # use first 10 seq if seq number larger than 10
            print(f"More than 10 sequences in full dataset. Automatically use first 10 sequences...\n")
            with open(subset_path + "_10", 'w') as file:
                seq = seq[0:10]
                seq_write = [f">{item}" for item in seq]  # add">" to each sequence
                file.write("".join(seq_write))
            print(f"Generating alignment plots for full dataset...\n")
            prettyplot(subset_path + "_10", subset)
        elif len(seq) <= 10 and len(seq) >= 0:
            print(f"Generating alignment plots for full dataset...\n")
            prettyplot(subset_path, subset)
print("All missions complete.\n")

