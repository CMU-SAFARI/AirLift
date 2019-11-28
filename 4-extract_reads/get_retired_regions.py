# python get_retired_regions.py [anticipated_read_size] [merged_file] [constant_regions] [output_file]
# anticipated_read_size is the size of the reads that you are mapping. 
# merged_regions: takes in a bed file of merged regions (from the regions that the gap sequences map to in the old reference)
# constant_regions: takes in chain file to identify constant regions. 
#
# output_file: the filename that dumps the bed file for retired regions. 

import sys 
import os.path 
#/home/firtinac/panzer/genome_remap/data/alignment_files/aligned_reads_bwa/merged/merged.bed

#https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def interval_merge(intervals): 
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

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
            else:
                merged.append(higher)
    return merged


chr_sizes = dict() 

    
    
anticipated_read_size = int(sys.argv[1]) # the size of the reads that you wil be mapping 
merged_regions = sys.argv[2]    # takes the merged bed file for all the reads that map to regions in the gap sequences. 
constant_regions = sys.argv[3]  # Takes a chain file
output_file = sys.argv[4]  # the file that is written to containing the retired regions in bed format. 


merged = dict() 
for line in open(merged_regions):
    chr1 = line.split()[0] 
    start = int(line.split()[1])
    end = int(line.split()[2]) 
    if chr1 not in merged:
        merged[chr1] = []
    merged[chr1].append([start, end]) 

constant = dict() 
chain = dict() 
chain_metadata = dict() 
print("Analyzing File: ", constant_regions) 
score = 0; tName = 0; tSize = 0; tStrand = 0; tStart = 0; tEnd = 0; qName = 0; qSize = 0; qStrand = 0; qStart = 0; qEnd = 0; c_id = 0; size = 0; dt = 0; dq = 0 
tstart_int = 0 

for line in open(constant_regions): 
    if line[0] == "#": 
        continue 
    if "chain" in line: 
        words = line.split() 
        score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, c_id = words[1:13] 
        tstart_int = int(tStart) 
        if tName not in chr_sizes: 
            chr_sizes[tName] = int(tSize)  
        if tName not in merged: 
            merged[tName] = [] 
    elif line in ['\n', '\r\n']:
        # do processing here 
        x = 1 # NOP 
    else: # size dt dq layout 
        words = line.split() 
        #print words 
        if len(words) == 3: 
            # calculate intervals for the data.  
            tend_int = tstart_int + (int(words[0]) - 1) # inclusive. 
            merged[tName].append([tstart_int, tend_int]) 
            tstart_int = tstart_int + int(words[0]) + int(words[1]) 
        else: 
            # calculate intervals for the data. 
            tend_int = tstart_int + (int(words[0]) - 1) 
            merged[tName].append([tstart_int, tend_int]) 

for qName in merged: 
    chr_size = 0 
    chr_size = chr_sizes[qName]
    
    for element in merged[qName]: 
        if element[0] >= anticipated_read_size: 
            element[0] = element[0] - anticipated_read_size
        else: 
            element[0] = 0 
        if element[1] <= chr_size - 1 - anticipated_read_size: 
            element[1] = element[1] + anticipated_read_size 
        else: 
            element[1] = element[1] = chr_size - 1 
    merged[qName] = interval_merge(merged[qName]) 

# Invert to get the retired regions. 
retired = dict() 
for qName in merged: 
    print("qName: ", qName) 
    size = 0 
    last_start = 0
    last_end = 0 
    elem_num = 0 
    if qName not in retired: 
        retired[qName] = [] 
    for element in merged[qName]: 

        if last_end > element[0]: 
            print("PROBLEM: " + str(last_end) + " " + element[0]) 
        size = size + element[1] - element[0] 

        if elem_num == 0:
            if element[0] > 0:  
                retired[qName].append([0, element[0]-1]) 
        else: 
            if (element[0] - 1) >= last_end+1: 
                retired[qName].append([last_end+1, element[0]-1])

        last_start = element[0]
        last_end = element[1] 
        elem_num = elem_num + 1 

    print("Size: ", size) 
    last_hole = -1 
    if len(retired[qName]) > 0: 
        last_hole = retired[qName][-1][1] 
    last_block = -1 
    if len(merged[qName]) > 0: 
        last_block = merged[qName][-1][1] 
    chr_size = chr_sizes[qName] 

    if len(retired[qName]) > 0: 
        if retired[qName][-1][1] < chr_size - 1: 
            if last_block != -1 and last_block < chr_size - 1: 
                retired[qName].append([last_block + 1, chr_size - 1]) 


ofn = open(output_file, "w") 
total = 0 
lines = 0 

for qName in merged: 
    for el in merged[qName]: 
        total = total + el[1] - el[0] 
        lines = lines + 1 
        #print qName + "\t" + str(el[0]) + "\t" + str(el[1]) 
print("merged -- total size: ", total) 
print("       -- lines     : ", lines) 

lines = 0 
total = 0 
for qName in retired: 
    for el in retired[qName]: 
        total = total + el[1] - el[0] 
        #print qName + "\t" + str(el[0]) + "\t" + str(el[1]) 
        write_string = qName + "\t" + str(el[0]) + "\t" + str(el[1]) + "\n" 
        ofn.write(write_string) 
        lines = lines + 1 
        
print("retired -- total size: ", total)
print("        -- lines     : ", lines) 
# TODO: merge gaps and merged. 
