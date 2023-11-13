#!/usr/bin/python3
import sys,subprocess,os
print("Starting plotcon...")
#winsize=input("Please enter the windowsize value(default:4): ") or 4

os.makedirs("./conservation_plot", exist_ok=True)
def plotco(winsize,path,i) :
  if i == 1:
     title="Full"
  elif i==2:
     title="Subset"
  plotcon_commands=[
                      f"plotcon {path} -winsize {winsize} -graph x11 -gsubtitle {title}", #I cannot prompt and save them at the same time.
                      f"plotcon {path} -winsize {winsize} -graph png -gdirectory ./conservation_plot -goutfile plot{i} "
                   ]
  for plotcon_command in plotcon_commands :
     plotcon=subprocess.Popen(plotcon_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     stdout, stderr = plotcon.communicate()
     print(stdout.decode())  # 打印标准输出结果
     print(stderr.decode())  # 打印错误输出结果
  return 0

fullset="./Clusout.fa"
subset="./subClusout.fa"
choice=input("Which dataset to use?(full/subset/both):") or "both"
winsize=input("Please enter the windowsize value(default:4): ") or 4
if choice == "both" :
  plotco(winsize,fullset,1)
  plotco(winsize,subset,2)
elif choice == "full" :
  plotco(winsize,fullset,1)
elif choice == "subset" :
  plotco(winsize,subset,2)
