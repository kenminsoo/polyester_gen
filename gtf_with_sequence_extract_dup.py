#gtf files may sometimes contain a sequence entry
#this will extract that entry along with the transcript treating the gtf file as a tsv
#also this is definitely easier done with bedtools, especially in the general case where sequences are not availabel
#which they are typically not

import pandas as pd
import subprocess
import random

class my_dictionary(dict):
 
  # __init__ function
  def __init__(self):
    self = dict()
 
  # Function to add key:value
  def add(self, key, value):
    self[key] = value

#set up some variables
out_dir = "out/"

subprocess.run(["mkdir", out_dir])

fasta_out = out_dir + "unique_transcripts.fasta"
summary_out = out_dir + "summary.txt"
gtf_path = "/Users/kenminsoo/Desktop/human_allRNA.gtf"
dup_remove = True
summary_sub_out = out_dir + "summary_sub.txt"
sub_gtf = out_dir + "sub_human_gtf.gtf"

#read in gtf
gtf_data = pd.read_table(gtf_path, header = None)
gtf_data.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf_data[['attribute', 'transcript_copy_id', 'sequence']] = gtf_data['attribute'].str.split(';',expand = True)

#keep only non-duplicated entries for base testings (i.e. the pipelines can function with no external variables)
all_entry = len(gtf_data)
#gtf_unique = gtf_data.drop_duplicates(subset=['attribute'], keep = False)
#gtf_unique = gtf_unique.drop_duplicates(subset=['sequence'], keep = False)
gtf_unique = gtf_data
num_unique = len(gtf_unique)

print(str(all_entry - num_unique) + " Transcripts have been removed")

#fast conversion
unique_fasta = gtf_unique[['attribute', 'sequence', 'source']]

#
num_5S = 0
num_rRNA = 0
num_transcripts = 0
databases = my_dictionary()

with open(fasta_out, 'w') as f:
    for index, row in unique_fasta.iterrows():
        
        transcript_name = ">" + str(row['attribute'])[15:-1] + "_" + str(row['source'])

        if "5S" in str(row['attribute'])[15:-1]:
            num_5S = num_5S + 1
        #    continue
        if "rRNA" in str(row['attribute'])[15:-1]:
            num_rRNA = num_rRNA + 1
        #    continue


        print(">" + str(row['attribute'])[15:-1] + "._." + str(row['source']), file=f)
        print(str(row['sequence'])[11:-1], file=f)

        if str(row['source']) not in databases:
            databases.add(str(row['source']), 0)

        databases[str(row['source'])] = databases[str(row['source'])] + 1

        num_transcripts = num_transcripts + 1

with open(summary_out, 'w') as f:
    print("number of 5S RNA removed: " + str(num_5S), file=f)
    print("number of rRNA removed: " + str(num_rRNA), file=f)
    print("number of transcripts in fasta: " + str(num_transcripts), file=f)
    print("number of RNA types in fasta: " + str(databases), file=f)

#we will take a random subset of a fasta file & create gtf file for the subset

subset_num = 1400

random_set = []

#set path
transcript_path = 'out/unique_transcripts.fasta'
subset_path = 'out/subset_transcripts.fasta'

transcript_list = []

#pick number of each transcripts based on the databases
#im not too sure if this is how i should do it...
#however, I think that this is what I need to do for now--taking user input
#may fix in future versions to not rely on human input as will be a pain for users
#but will be integrating more datasets with the annotation file so i'd like to keep this
#at least some what applicable if things begin to change 
set_numb_types = my_dictionary()

for set in databases:
    print(set)
    num_of_trans = input("How many transcripts would you like for the printed set:")
    num_of_trans = int(num_of_trans)
    set_numb_types.add(set, num_of_trans)     

#i can't tell if this a really bad way to randomly pick transcripts of x type
#it should work at an okay speed though, it is doubtful that users will pick a subset of 40k
source_sum = 0
with open(transcript_path, 'r') as f:
    contents = f.readlines()
#here, we populate the random number list with each transcript type, one by one, till number is reached
    for source in databases:
        source_sum = source_sum + set_numb_types[source]
        print(source_sum)
        while len(random_set) < source_sum:
            the_number = random.randint(0, num_transcripts-1)
            the_number = (the_number * 2)
            if source in contents[the_number].strip():
                if the_number not in random_set:
                    random_set.append(the_number)

    for num in random_set:
        transcript_list.append(contents[num].strip())
        transcript_list.append(contents[num+1].strip())

subset_summary = my_dictionary()
rRNA_types = {"5S":0, "rRNA":0}

with open(subset_path, 'w')  as f:
    for item in transcript_list:
        for source in databases:
            if source in item: 
                if source not in subset_summary:
                    subset_summary.add(source, 0)
                subset_summary[source] = subset_summary[source] + 1
        #count the rRNA types
        if "rRNA" in item:
            if "rRNA" not in rRNA_types:
                rRNA_types = rRNA_types.add("rRNA", 0)
            rRNA_types["rRNA"] = rRNA_types["rRNA"] + 1
        if "5S" in item:
            if "5S" not in rRNA_types:
                rRNA_types = rRNA_types.add("5S", 0)
            rRNA_types["5S"] = rRNA_types["5S"] + 1
        print(item.partition('._.')[0], file=f)

#with open(summary_sub_out, 'w') as f:
    #print(subset_summary, file=f)
    #print(rRNA_types, file=f)

#note that the following gtf file creation will only work if each transcript is unique with unique sequence
#length_fast = len(transcript_list)
#for i in range(int(length_fast)):
    #if i%2 == 1:
        #transcript_list[i]
    #else:
        #continue