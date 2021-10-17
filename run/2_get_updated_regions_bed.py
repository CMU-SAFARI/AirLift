# extract the updated bed regions from the exact chain file 
# these regions will be used to construct the CIGAR chain file 
# Assumes that no negative strands occur on the old references 


import sys
import os.path 

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
            if higher[0] < lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
    return merged

if not len(sys.argv) == 5: 
    print("Not enough arguments.") 
    print("Usage: python analyze_liftover.py <chain_file> <read size> <max errors in alignment> <output bedfile>") 
    exit() 

if not os.path.isfile(sys.argv[1]): 
    exit()

verbose = 0

fn_chain = sys.argv[1] 
read_size = int(sys.argv[2])
max_errors = int(sys.argv[3]) 
output_file = sys.argv[4]

chain = dict() 

if (verbose): print("Analyzing File: ", fn_chain)

score = 0; tName = 0; tSize = 0; tStrand = 0; tStart = 0; tEnd = 0; qName = 0; qSize = 0; qStrand = 0; qStart = 0; qEnd = 0; c_id = 0; size = 0; dt = 0; dq = 0 
chain_region_chr = 0
chain_region_start = 0
chain_region_end = 0  
gap_start = 0
gap_end = 0 

chr_sizes = dict() 

# immediately build gaps in chain dict 
# gaps are inclusive [] 
for line in open(fn_chain): 
    fields = line.split() 
    if len(fields) > 0 and fields[0] == "#": 
        continue 

    if "chain" in fields: 
        score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, c_id = fields[1:13] 

        if tName not in chr_sizes: 
            chr_sizes[tName] = int(tSize) 
        if qName not in chr_sizes:
            chr_sizes[qName] = int(qSize) 
        chain_region_chr = tName 
        chain_region_start = int(tStart) 
        chain_region_end = int(tEnd) 
        gap_start = chain_region_start 
        gap_end = chain_region_start 

        if chain_region_chr not in chain:        
            chain[chain_region_chr] = [] 

    elif len(fields) == 3: 
        gap_start = gap_start + int(fields[0]) 
        gap_end = gap_start   + int(fields[1]) - 1 

        # found the gap 
        chain[chain_region_chr].append([gap_start, gap_end]) 
        
        # reset gap start 
        gap_start = gap_start + int(fields[1]) 
        
    elif len(fields) == 1:  
        # do nothing. 
        x=0 

if (verbose): print(chain) 

# extend the gaps 
for name in chain: 
    chr_size = chr_sizes[name] 

    for el in chain[name]: 
        if el[0] >= (read_size + max_errors): 
            el[0] = el[0] - (read_size + max_errors)
        else: 
            el[0] = 0 
        if el[1] <= chr_size - 1 - (read_size + max_errors): 
            el[1] = el[1] + (read_size + max_errors) 
        else: 
            el[1] = chr_size - 1 

if (verbose): print(chain) 

# merge the gaps 
for name in chain: 
    chain[name] = interval_merge(chain[name]) 
        
if (verbose): print(chain) 


# write the intervals to a file 
f = open(output_file, "w") 
for name in chain: 
    for el in chain[name]: 
        f.write(name + "\t" + str(el[0]) + "\t" + str(el[ 1]) + "\n") 
f.close() 
    
