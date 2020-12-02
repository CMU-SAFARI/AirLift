# python get_constant_regions.py [anticipated_read_size] [chain file] [output bed file]
# anticipated_read_size is the size of the reads that you are mapping. 
# chain file ths chain file between old and new reference genome 
# OUTPUT:  bed file of the constant regions between the reference genomes

# How to Run Example: python get_constant_regions.py 101 /mnt/panzer/firtinac/genome_remap/data/chain_files/hg19tohg38/hg19tohg38_complete.chain constant.bed 

# TODO: change the filename of merged file to include the hgXtohgX 

# when using a liftover file Hg19to38. 19 is the target, and 38 is the query reference genome. 

import sys
#import numpy as np 
import os.path 
import math 
import bisect 

#https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def interval_merge(intervals): 
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    overlap = 0 

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
                if higher[1] <= lower[1]: 
                    overlap = overlap + higher[1] - higher[0] + 1 
                else: 
                    overlap = overlap + lower[1] - higher[0] + 1 
            else:
                merged.append(higher)
    return merged, overlap 

old_chr_sizes = dict() 
new_chr_sizes = dict() 

anticipated_read_size = int(sys.argv[1]) # the size of the reads that you wil be mapping 
ofn = open(sys.argv[3], 'w')  

files = dict() 
merged = dict()
constant = dict()
chain = dict() 
chain_metadata = dict() 
#alignment_file = sys.argv[2] 
score = 0; tName = 0; tSize = 0; tStrand = 0; tStart = 0; tEnd = 0; qName = 0; qSize = 0; qStrand = 0; qStart = 0; qEnd = 0; c_id = 0; size = 0; dt = 0; dq = 0 

args = sys.argv[2] 
fn = args.split("/")[-1] 
# parse as a chain file 
if ".chain" in fn: 
    references = fn.split("_")[0] 
    if references not in files: 
        files[references] = ["","",""] 
    files[references][2] = args 

# parse as a bed file 
elif ".bed" in fn: 
    references = fn.split("_")[0] 
    if references not in files: 
        files[references] = ["","",""] 
    files[references][0] = args 

elif os.path.isdir(args): 
    references = args.split("/")[-1].split("_")[-1] 
    if references not in files: 
        files[references] = ["","",""] 
    files[references][1] = args 

else: 
    print(args + " could not be parsed") 

#print files 

a = dict() 

for references in files: 
    tstart_int = 0
    qstart_int = 0 
    tend_int = 0
    qend_int = 0 
    old_chr_sizes = dict() 
    new_chr_sizes = dict() 
    print(references) 
    a[references] = [0,0,0,0] 
    chain = dict() 
    chain_metadata = dict() 
    merged = dict() 
    retired = dict() 
    new = dict() 
    real_new = dict() 
    constant = dict() 
    updated = dict() 
    constant_region_size = 0 
    updated_size = 0 
    overlap_old_constant_updated = 0 
    overlap_new_constant_updated = 0 
    for fn in files[references]: 
        if ".chain" in fn: 
            #print fn 
            for line in open(fn): 
                if line[0] == "#": 
                    continue 
                if "chain" in line: 
                    words = line.split() 
                    score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, c_id = words[1:13] 
                    qstart_int = int(qStart) 
                    tstart_int = int(tStart) 
                    if tName not in constant: 
                        constant[tName] = [] 
                elif line in ['\n', '\r\n']:
                    # do processing here 
                    x = 1 # NOP 
                else: # size dt dq layout 
                    words = line.split() 
                    #print words 
                    if len(words) == 3: 
                        # calculate intervals for the data.  
                        tend_int = tstart_int + (int(words[0]) - 1) # inclusive. 
                        qend_int = qstart_int + (int(words[0]) - 1) 
                        constant[tName].append([tstart_int, tend_int]) 
                        constant_region_size = constant_region_size + tend_int - tstart_int 
                        tstart_int = tstart_int + int(words[0]) + int(words[1]) 
                        qstart_int = qstart_int + int(words[0]) + int(words[2]) 
                    else: 
                        # calculate intervals for the data. 
                        tend_int = tstart_int + (int(words[0]) - 1) 
                        qend_int = qstart_int + (int(words[0]) - 1) 
                        constant[tName].append([tstart_int, tend_int]) 
                        constant_region_size = constant_region_size + tend_int - tstart_int 

lines = 0 
total = 0 
for qName in constant: 
    for el in constant[qName]: 
        total = total + el[1] - el[0] 
        #print qName + "\t" + str(el[0]) + "\t" + str(el[1]) 
        write_string = qName + "\t" + str(el[0]) + "\t" + str(el[1]) + "\n" 
        ofn.write(write_string) 
        lines = lines + 1 
