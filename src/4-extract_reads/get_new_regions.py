# python get_new_regions.py [anticipated_read_size] [merged_file] [constant_regions] [output_file] 
# anticipated_read_size is the size of the reads that you are mapping. 
# merged_regions: takes in a bed file of merged regions (from the regions that the gap sequences map to in the old reference)
# constant_regions: takes in chain file to identify constant regions. 

# How to Run Example:  python region_composition_analysis.py 101 /mnt/panzer/firtinac/genome_remap/data/chain_files/hg19tohg38/hg19tohg38_complete.chain /mnt/panzer/firtinac/genome_remap/data/alignment_files/aligned_reads_bwa_hg19tohg38/merged/hg19tohg38_merged.bed /mnt/panzer/firtinac/genome_remap/data/alignment_files/aligned_reads_bwa_hg19tohg38

# TODO: change the filename of merged file to include the hgXtohgX 

# when using a liftover file Hg19to38. 19 is the target, and 38 is the query reference genome. 

import sys
import matplotlib.pyplot as plt
import numpy as np 
import os.path 
import math 

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

files = dict() 
merged = dict()
constant = dict()
chain = dict() 
chain_metadata = dict() 
score = 0; tName = 0; tSize = 0; tStrand = 0; tStart = 0; tEnd = 0; qName = 0; qSize = 0; qStrand = 0; qStart = 0; qEnd = 0; c_id = 0; size = 0; dt = 0; dq = 0 

for args in sys.argv[2:]: 
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
        print args + " could not be parsed" 

#print files 

a = dict() 

for references in files: 
    tstart_int = 0
    qstart_int = 0 
    tend_int = 0
    qend_int = 0 
    old_chr_sizes = dict() 
    new_chr_sizes = dict() 
    print references 
    a[references] = [0,0,0,0] 
    chain = dict() 
    chain_metadata = dict() 
    merged = dict() 
    retired = dict() 
    new = dict() 
    real_new = dict() 
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
                    if tName not in old_chr_sizes: 
                        old_chr_sizes[tName] = int(tSize)  
                    if qName not in new_chr_sizes: 
                        new_chr_sizes[qName] = int(qSize) 
                    if tName not in merged: 
                        merged[tName] = [] 
                    if qName not in new: 
                        new[qName] = [] 
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
                        merged[tName].append([tstart_int, tend_int]) 
                        new[qName].append([qstart_int, qend_int])
                        constant_region_size = constant_region_size + tend_int - tstart_int 
                        tstart_int = tstart_int + int(words[0]) + int(words[1]) 
                        qstart_int = qstart_int + int(words[0]) + int(words[2]) 
                    else: 
                        # calculate intervals for the data. 
                        tend_int = tstart_int + (int(words[0]) - 1) 
                        qend_int = qstart_int + (int(words[0]) - 1) 
                        merged[tName].append([tstart_int, tend_int]) 
                        new[qName].append([qstart_int, qend_int])
                        constant_region_size = constant_region_size + tend_int - tstart_int 

            for qName in merged: 
                chr_size = 0 
                chr_size = old_chr_sizes[qName]

                for element in merged[qName]: 
                    if element[0] >= anticipated_read_size: 
                        element[0] = element[0] - anticipated_read_size
                    else: 
                        element[0] = 0 
                    if element[1] <= chr_size - 1 - anticipated_read_size: 
                        element[1] = element[1] + anticipated_read_size 
                    else: 
                        element[1] = element[1] = chr_size - 1 
                merged[qName], x = interval_merge(merged[qName]) 
                overlap_new_constant_updated = overlap_new_constant_updated + x 
            
            for qName in new: 
                chr_size = 0 
                chr_size = new_chr_sizes[qName]

                for element in new[qName]: 
                    if element[0] >= anticipated_read_size: 
                        element[0] = element[0] - anticipated_read_size
                    else: 
                        element[0] = 0 
                    if element[1] <= chr_size - 1 - anticipated_read_size: 
                        element[1] = element[1] + anticipated_read_size 
                    else: 
                        element[1] = element[1] = chr_size - 1 
                new[qName], x = interval_merge(new[qName]) 
                overlap_old_constant_updated = overlap_old_constant_updated + x 

# Invert to get the real new regions. 
            for qName in new: 
                #print "qName: ", qName 
                size = 0 
                last_start = 0
                last_end = 0 
                elem_num = 0 
                if qName not in real_new: 
                    real_new[qName] = [] 
                for element in new[qName]: 

                    if last_end > element[0]: 
                        print "PROBLEM: ", last_end, " ", element[0] 
                    size = size + element[1] - element[0] 

                    if elem_num == 0:
                        if element[0] > 0:  
                            real_new[qName].append([0, element[0]-1]) 
                    else: 
                        if (element[0] - 1) >= last_end+1: 
                            real_new[qName].append([last_end+1, element[0]-1])

                    last_start = element[0]
                    last_end = element[1] 
                    elem_num = elem_num + 1 

                #print "Size: ", size 
                last_hole = -1 
                if len(real_new[qName]) > 0: 
                    last_hole = real_new[qName][-1][1] 
                last_block = -1 
                if len(new[qName]) > 0: 
                    last_block = new[qName][-1][1] 
                chr_size = old_chr_sizes[qName] 

                if len(real_new[qName]) > 0: 
                    if real_new[qName][-1][1] < chr_size - 1: 
                        if last_block != -1 and last_block < chr_size - 1: 
                            real_new[qName].append([last_block + 1, chr_size - 1]) 

# Invert to get the retired regions. 
            for qName in merged: 
                #print "qName: ", qName 
                size = 0 
                last_start = 0
                last_end = 0 
                elem_num = 0 
                if qName not in retired: 
                    retired[qName] = [] 
                for element in merged[qName]: 

                    if last_end > element[0]: 
                        print "PROBLEM: ", last_end, " ", element[0] 
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

                #print "Size: ", size 
                last_hole = -1 
                if len(retired[qName]) > 0: 
                    last_hole = retired[qName][-1][1] 
                last_block = -1 
                if len(merged[qName]) > 0: 
                    last_block = merged[qName][-1][1] 
                chr_size = old_chr_sizes[qName] 

                if len(retired[qName]) > 0: 
                    if retired[qName][-1][1] < chr_size - 1: 
                        if last_block != -1 and last_block < chr_size - 1: 
                            retired[qName].append([last_block + 1, chr_size - 1]) 

        elif ".bed" in fn: 
            updated_size = 0 
            for line in open(fn):
                chr1 = line.split()[0] 
                start = int(line.split()[1])
                end = int(line.split()[2]) 
                if chr1 not in merged:
                    merged[chr1] = []
                merged[chr1].append([start, end]) 
                updated_size = updated_size + end - start 
        
        elif os.path.isdir(fn):
            # read in all the .bed file names in this path and merge it with the news. 
            for filename in os.listdir(fn):
                if ".bed" in filename: 
                    chrname = filename.split("_")[0] 
                    start = int(filename.split("_")[1])
                    end = int(filename.split("_")[2])
                    if chrname not in new: 
                        new[chrname] = [] 
                    new[chrname].append([start, end])

ofn = open(output_file, "w") 

total_size = 0 
for qName in real_new:  
    for el in real_new[qName]: 
        total_size = total_size + el[1] - el[0] 
        write_string = qName + "\t" + str(el[0]) + "\t" + str(el[1]) + "\n" 
        ofn.write(write_string) 

print "New Region total_size: ", total_size 
    

#    lines = 0 
#    retired_size = 0 
#    for qName in retired: 
#        for el in retired[qName]: 
#            retired_size = retired_size + el[1] - el[0] 
#            #print qName + "\t" + str(el[0]) + "\t" + str(el[1]) 
#            #write_string = qName + "\t" + str(el[0]) + "\t" + str(el[1]) + "\n" 
#            #ofn.write(write_string) 
#            lines = lines + 1 
#
#    new_size = 0 
#    for qName in real_new: 
#        for el in real_new[qName]: 
#            new_size = new_size + el[1] - el[0] 
#
#    old_chr_size = 0 
#    for c in old_chr_sizes: 
#        old_chr_size = old_chr_size + old_chr_sizes[c] 
#        
#    new_chr_size = 0 
#    for c in new_chr_sizes:
#        new_chr_size = new_chr_size + new_chr_sizes[c] 
#
#    print "References:   ", references 
#    print "old chr size: ", str(old_chr_size)  
#    print "new chr size: ", str(new_chr_size)  
#    print "between new constant and updated overlap: ", overlap_new_constant_updated 
#    print "between old constant and updated overlap: ", overlap_old_constant_updated 
#    
#    print "updated  -- total size: ", updated_size 
#    print "retired  -- total size: ", retired_size
##    print "         -- lines     :  ", lines 
#         
#    print "constant -- total size: ", constant_region_size 
#    print "new      -- total size: ", new_size 
#    a[references] = [constant_region_size, updated_size, retired_size, new_size] 
#
##print a 
#for ref in a: 
#    print ref
#    total_size = float(sum(a[ref])) 
#    print "   - constant size: " + str(a[ref][0]) + " (" + str(math.ceil(a[ref][0]/total_size*1000)/1000) + ")" 
#    print "   -  updated size: " + str(a[ref][1]) + " (" + str(math.ceil(a[ref][1]/total_size*1000)/1000) + ")" 
#    print "   -  retired size: " + str(a[ref][2]) + " (" + str(math.ceil(a[ref][2]/total_size*1000)/1000) + ")" 
#    print "   -      new size: " + str(a[ref][3]) + " (" + str(math.ceil(a[ref][3]/total_size*1000)/1000) + ")" 
#
#print "====================== REFORMATTED ========================" 
#
#for ref in a: 
#    print ref 
#    total_size = float(sum(a[ref])) 
#    print str(math.ceil(a[ref][0]/total_size*1000)/1000) + "/" + str(math.ceil(a[ref][1]/total_size*1000)/1000) + "/" + str(math.ceil(a[ref][2]/total_size*1000)/1000) + "/" + str(math.ceil(a[ref][3]/total_size*1000)/1000)
#
##ax, fig = plt.subplots() 



