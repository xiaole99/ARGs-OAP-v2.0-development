
#!/usr/bin/env python
from Bio import SeqIO 
import os
import sys
import csv

filename=raw_input('Pls enter the full name of type.list: ')
filename2=raw_input('Pls enter the full name of extracted_sequence_fasta.file: ')
foldername=raw_input('Pls enter the pathway of folder: ')

fileoutput1=open("single_keyword_list.txt","w")
fileoutput2=open("keywords_with_space.txt","w")

f1=open(filename,'r')
i=0
j=0
for line in f1:
    if " " in line:
        i+=1
        fileoutput2.write(str(line).strip()+"\n")
    else:
        j+=1
        fileoutput1.write(str(line).strip()+"\n")

fileoutput1.close()
fileoutput2.close()

print j, "single keywords are listed!"
print i, "keywords has space inside!"

a={}
for record in SeqIO.parse(filename2,"fasta"):
    a[">"+str(record.description).strip()+"\n"+str(record.seq)]=1
A={}
for record in SeqIO.parse(filename2,"fasta"):
    A[str(record.id).strip()]=1
    
file1=open("single_keyword_list.txt","r")

b={}
NUM=0
subtype={}
SUBTYPE={}
for line in file1:
    Num=0
    subtype[str(line).strip()]=[]
    SUBTYPE[str(line).strip()]=[]
    for record in SeqIO.parse(filename2,"fasta"):
        Num +=1
        if str(line).strip().lower() in str(record.description).lower().split():
            if ">"+str(record.description).strip()+"\n"+str(record.seq) in a.keys():
                subtype[str(line).strip()].append(">"+str(record.description).strip()+"\n"+str(record.seq))
                SUBTYPE[str(line).strip()].append(str(record.id).strip())
                del a[">"+str(record.description).strip()+"\n"+str(record.seq)]
                del A[str(record.id).strip()]
            else:
                b[">"+str(record.description).strip()+"\n"+str(record.seq)]= 1
                for key,value in subtype.items():
                    if ">"+str(record.description).strip()+"\n"+str(record.seq) in value:
                        subtype[key].remove(">"+str(record.description).strip()+"\n"+str(record.seq))
                for key,value in SUBTYPE.items():
                    if str(record.id).strip() in value:
                        SUBTYPE[key].remove(str(record.id).strip())
    print Num, "sequences have been searched for_{}".format(str(line).strip()), "to match single keywords"
 

file2=open("keywords_with_space.txt","r")
n=0
for line in file2:
    N=0
    subtype[str(line).strip()]=[]
    SUBTYPE[str(line).strip()]=[]
    for record in SeqIO.parse(filename2,"fasta"):
        N +=1
        if str(line).strip().lower() in str(record.description).lower():
            if ">"+str(record.description).strip()+"\n"+str(record.seq) in a.keys():
                subtype[str(line).strip()].append(">"+str(record.description).strip()+"\n"+str(record.seq))
                SUBTYPE[str(line).strip()].append(str(record.id).strip())
                del a[">"+str(record.description).strip()+"\n"+str(record.seq)]
                del A[str(record.id).strip()]
            else:
                b[">"+str(record.description).strip()+"\n"+str(record.seq)]= 1
                for key,value in subtype.items():
                    if ">"+str(record.description).strip()+"\n"+str(record.seq) in value:
                        subtype[key].remove(">"+str(record.description).strip()+"\n"+str(record.seq))
                for key,value in SUBTYPE.items():
                    if str(record.id).strip() in value:
                        SUBTYPE[key].remove(str(record.id).strip())
    print N, "sequences have been searched for_{}".format(str(line).strip()), "to match keywords with space"
print len(b), "replicates are there!"

print len(a), "sequences have not been matched!"
print len(A), "IDs have not been matched!"

fileoutput3=open("total_"+str(len(a))+"_notmatch.txt", "w") 
fileoutput4=open("total_"+str(len(b))+"_replicate.txt","w")
for key in a.keys():
    fileoutput3.write(str(key)+"\n")
fileoutput3.close()

for key in b.keys():
    fileoutput4.write(str(key)+"\n")
fileoutput4.close()
   

for key in subtype.keys():
    if "/" in key:
        subtype[key.replace("/","_")]= subtype.pop(key)
        
for k, lists in subtype.items():
    with open(os.path.join(foldername,"{}.txt".format(k)), "w") as out:
        out.write('\n'.join(lists)+'\n')

with open('structure20161102.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file,delimiter='\t')
    for key, value in SUBTYPE.items():
       writer.writerow([key,value])
        

print 'OK, Finished!'
raw_input('Press <Enter> to close this window!')


