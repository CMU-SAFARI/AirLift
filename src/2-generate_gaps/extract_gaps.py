# python extract_gaps.py <chain file> <reference genome FASTA> <anticipated_read_size> <output FASTA file of gaps> 
# just look at the liftover file <arg1> and list all the locations of gaps across the chain file and the length of the gaps. 
# the reference genome FASTA file <arg2> that the liftover file targets (Query in the liftover chain file)  
# 50 base pairs per line? 
# anticipated_read_size is the size of the reads that you want to map to these gaps. This just extends each gap both sides by anticipated_read_size and prints out each of these gaps so we can map to just them 
# outputs a file in the form of FASTA of all the gaps in the chain file, so we can map the reads to this FASTA. 

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
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
    return merged

# stats 
num_chains = 0 

if not len(sys.argv) == 5: 
    print("Not enough arguments.") 
    print("Usage: python analyze_liftover.py <chain_file> <reference genome FASTA file> <anticipated read size> <output FASTA FILE of gaps>") 
    exit() 

if not os.path.isfile(sys.argv[1]): 
    exit()

fn_chain = sys.argv[1] 
anticipated_read_size = int(sys.argv[3])

if os.path.isfile(sys.argv[4]): 
    prompt = "overwrite " + sys.argv[4] + "? [Enter] for yes [any other key] for no" 
    text = input(prompt) 
    if text != "":
        exit() 

chain = dict() 
chain_metadata = dict() 
print("Analyzing File: ", fn_chain)
score = 0; tName = 0; tSize = 0; tStrand = 0; tStart = 0; tEnd = 0; qName = 0; qSize = 0; qStrand = 0; qStart = 0; qEnd = 0; c_id = 0; size = 0; dt = 0; dq = 0 
tstart_int = 0 
tend_int = 0 

for line in open(fn_chain): 
    if line[0] == "#":
        continue 
    if "chain" in line: 
        words = line.split() 
        score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, c_id = words[1:13] 
        tstart_int = int(tStart) 
        qstart_int = int(qStart) 
        qChain_offset = 0 
        tChain_offset = 0 
        if qName not in chain: 
            chain[qName] = dict() 
        if c_id not in chain[qName]: 
            chain[qName][c_id] = dict() 
            chain[qName][c_id]["meta"] = [int(score), tName, int(tSize), tStrand, int(tStart), int(tEnd), qName, int(qSize), qStrand, int(qStart), int(qEnd), int(c_id)] 
            chain_metadata[qName] = dict() 
            chain_metadata[qName][c_id] = [int(qEnd) - int(qStart),0,0,0] # size, gaps, gaps_size, blocks 
            chain[qName][c_id]["tdata"] = [] 
            chain[qName][c_id]["qdata"] = [] 
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
            chain[qName][c_id]["tdata"].append([tstart_int, tend_int]) 
            chain[qName][c_id]["qdata"].append([qstart_int, qend_int]) 
            tstart_int = tstart_int + int(words[0]) + int(words[1]) 
            qstart_int = qstart_int + int(words[0]) + int(words[2]) 

            chain_metadata[qName][c_id][1] = chain_metadata[qName][c_id][1] + 1 
            chain_metadata[qName][c_id][2] = chain_metadata[qName][c_id][2] + int(words[2]) 
            chain_metadata[qName][c_id][3] = chain_metadata[qName][c_id][3] + int(words[0]) 
        else: 
            # calculate intervals for the data. 
            tend_int = tstart_int + (int(words[0]) - 1) 
            qend_int = qstart_int + (int(words[0]) - 1) 
            chain[qName][c_id]["tdata"].append([tstart_int, tend_int]) 
            chain[qName][c_id]["qdata"].append([qstart_int, qend_int]) 

            chain_metadata[qName][c_id][3] = chain_metadata[qName][c_id][3] + int(words[0]) 

for qName in chain_metadata: 
    print(qName) 
    total_size = 0
    total_gaps = 0
    total_gapsize = 0 
    total_blocks = 0 
    for c_id in chain_metadata[qName]: 
        total_size = total_size + chain_metadata[qName][c_id][0] 
        total_gaps = total_gaps + chain_metadata[qName][c_id][1] 
        total_gapsize = total_gapsize + chain_metadata[qName][c_id][2] 
        total_blocks = total_blocks + chain_metadata[qName][c_id][3] 
    print( "total_size: ", total_size )
    print( "total_gaps: ", total_gaps )
    print( "total_gapsize: ", total_gapsize )
    print( "total_blocks: ", total_blocks ) 

# Logic block for checking overlaps across the query chromosome: There are. 
#for qName in chain: 
#    for c_id in chain[qName]: 
#        chromosome1 = chain[qName][c_id]["meta"][6] 
#        first_start = chain[qName][c_id]["meta"][9] # qStart 
#        first_end = chain[qName][c_id]["meta"][10] - 1 # qEnd non-inclusive. 
#        for c_id2 in chain[qName]: 
#            if c_id == c_id2:
#                continue 
#            chromosome2 = chain[qName][c_id]["meta"][6] 
#            if chromosome1 != chromosome2: # make sure we are only checking within chromosome 
#                continue 
#            second_start = chain[qName][c_id]["meta"][9] # qStart 
#            second_end = chain[qName][c_id]["meta"][10] - 1 # qEnd non-inclusive. 
#            if ((first_start >= second_start and first_start <= second_end) or 
#               (first_end >= second_start and first_end <= second_end)): 
#                print "overlap" 

# this is the list of non-overlapping gaps per chromosome. 
q_all_gaps = dict() # it starts off as the non-gaps, but then we invert it. 
gaps = dict() 


for qName in chain:
    print("chr: ", qName) 
    if qName not in q_all_gaps:    
        q_all_gaps[qName] = [] 
        gaps[qName] = [] 
    for c_id in chain[qName]: 
        for interval in chain[qName][c_id]["qdata"]: 
            q_all_gaps[qName].append(interval) 
    q_all_gaps[qName] = interval_merge(q_all_gaps[qName]) 

size = 0 
size_gaps = 0 

print( "=============") 
print( "BLOCKS: ") 
print( "=============") 

for qName in q_all_gaps: 
    print("qName: ", qName) 
    size = 0 
    last_start = 0
    last_end = 0 
    elem_num = 0 
    for element in q_all_gaps[qName]: 

        if last_end > element[0]: 
            print("PROBLEM: ", last_end, " ", element[0]) 
        size = size + element[1] - element[0] 

        if elem_num == 0:
            if element[0] > 0:  
                gaps[qName].append([0, element[0]-1]) 
        else: 
            if (element[0] - 1) >= last_end+1: 
                gaps[qName].append([last_end+1, element[0]-1])

        last_start = element[0]
        last_end = element[1] 
        elem_num = elem_num + 1 


    print("Size: ", size) 
    last_hole = -1 
    if len(gaps[qName]) > 0: 
        last_hole = gaps[qName][-1][1] 
    last_block = -1 
    if len(q_all_gaps[qName]) > 0: 
        last_block = q_all_gaps[qName][-1][1] 
    chr_size = 0 
    for c_id in chain[qName]: 
        chr_size = chain[qName][c_id]["meta"][7]
        break 

    if len(gaps[qName]) > 0: 
        if gaps[qName][-1][1] < chr_size - 1: 
            if last_block != -1 and last_block < chr_size - 1: 
                gaps[qName].append([last_block + 1, chr_size - 1]) 

print( "=============" )
print( "GAPS: ")
print( "=============" ) 

for qName in gaps: 
    print(qName) 
    size_gaps = 0 
    for element in gaps[qName]: 
       size_gaps = size_gaps + element[1] - element[0] + 1  

    print(size_gaps) 


# ===================== Do some plotting with our data. ======================
# first look at the CDF of the gaps across the entire chromosome. 
#for qName in gaps: 
#    x = [0] 
#    y = [0] 
#    if "_" in qName: # skip all alts and randoms. 
#        continue 
#    plt.figure() 
#    last_start = -1 
#    last_end = -1
#    for element in gaps[qName]: 
#        x.append(element[0])
#        y.append(y[-1]) 
#        x.append(element[1]) 
#        y.append(y[-1] + element[1]) 
#    
#    plt.plot(x,y) 
#    plt.title(qName) 
#
#plt.show() 

######################################################
# extend and merging gaps 
######################################################

num_gaps = 0 
for qName in gaps: 
    num_gaps = num_gaps + len(gaps[qName])

print("num_gaps before gap_extend_merge: ", num_gaps) 

# extend each gap by anticipated_read_size on either end and merge all gaps again. 
for qName in gaps: 
    chr_size = 0 
    for c_id in chain[qName]:
        chr_size = chain[qName][c_id]["meta"][7] 
        break 

    for element in gaps[qName]: 
        if element[0] >= anticipated_read_size: 
            element[0] = element[0] - anticipated_read_size
        else: 
            element[0] = 0 
        if element[1] <= chr_size - 1 - anticipated_read_size: 
            element[1] = element[1] + anticipated_read_size 
        else: 
            element[1] = chr_size - 1 

    gaps[qName] = interval_merge(gaps[qName]) 

num_gaps = 0 
for qName in gaps: 
    num_gaps = num_gaps + len(gaps[qName])
print("num_gaps after gap_extend_merge: ", num_gaps) 


########################################################
# splitting up the reference genome into the gaps 
########################################################

fasta_fn = open(sys.argv[2]) 
output_fn = open(sys.argv[4], "w") 

gap_offset = 0 
this_chr = "" 
bp_num = 0 
elem_num = 0
goto_next_chr = 0 
gap = [0,0] 
for line in fasta_fn: 
    if ">" in line: 
        this_chr = line.split(">")[1].split()[0] 
        print(this_chr) 
        bp_num = 0 
        elem_num = 0 
        goto_next_chr = 0 
        if this_chr not in gaps.keys(): 
            print(this_chr, " not in gaps..") 
            goto_next_chr = 1 
            continue 
        if elem_num < len(gaps[this_chr]): 
            gap = gaps[this_chr][elem_num] 
        
    else: 
        if goto_next_chr: 
            continue 
        for bp in line: 
            if bp.isspace(): 
                continue 

            gap_offset = bp_num - gap[0] 

            # 0 indexed. 
            if bp_num < gap[0]:
                bp_num = bp_num + 1 
                continue 

            if bp_num == gap[0]: 
                title = ">" + str(this_chr) + "_" + str(gap[0]) + "_" + str(gap[1]) + "_g" + str(elem_num) + "\n" 
                output_fn.write(title) 

            if bp_num >= gap[0] and bp_num <= gap[1]: 
                if (gap_offset > 0) and ((gap_offset) % 50 == 0) and (bp_num < gap[1]): 
                    output_fn.write("\n") 
                output_fn.write(bp) 

            if bp_num == gap[1]: 
                if len(gaps[this_chr]) - 1 > elem_num: 
                    elem_num = elem_num + 1 
                    gap = gaps[this_chr][elem_num] 
                else: 
                    goto_next_chr = 1 
                output_fn.write("\n") 
            
            bp_num = bp_num + 1 


