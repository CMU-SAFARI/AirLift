from Bio import SeqIO 
import sys 

# assumes that only the new regions will be negative stranded.. 

# running instructions: 
# to verify whether a chain file is an exact chain file: 
### python exact_chain_file.py verify old_ref.fa new_ref.fa chain_file 

# to build an exact chain file out of an inaccurate chain file 
### python exact_chain_file.py build old_ref.fa new_ref.fa chain_file 

# this indicates the smallest region size acceptable to be placed in a chain file. 
min_region_size = 100

# verbosity level. set to 1 to debug. 
verbose = 0 

rev_comp = {"c": "g", "g": "c", "a":"t", "t":"a", "n":"n"}

def compare_seqs(a, b): 
    if len(a) != len(b):
        return False 
    for el in range(len(a)): 
        if a[el].lower() != b[el].lower():
            return False 

    return True 

# B is rev comped 
def compare_revcomp_seqs(a, b): 
    if len(a) != len(b): 
        return False 
    blen = len(b) 
    for el in range(len(a)): 
        if a[el].lower() != rev_comp[b[blen - 1 - el].lower()]: 
            return False 
    return True


def read_chain_file(cf_name): 
    score = ""
    source_name = ""
    source_size = 0
    source_strand = ""
    source_start = 0
    source_end = 0 
    target_name = ""
    target_size = 0
    target_strand = ""
    target_start = 0
    target_end = 0
    sfrom = 0
    tfrom = 0 
     
    target_chromSize = {}
    source_chromSize = {} 
    maps = {} 
    chain_id = "" 

    for line in open(cf_name): 
        fields = line.split() 

        if "chain" in line: 
            score = int(fields[1])
            source_name = fields[2] 
            source_size = int(fields[3]) 
            source_strand = fields[4] 
            source_start = int(fields[5]) 
            source_end = int(fields[6]) 
            target_name = fields[7] 
            target_size = int(fields[8]) 
            target_strand = fields[9] 
            target_start = int(fields[10])
            target_end = int(fields[11]) 
            target_chromSize[target_name] = target_size
            source_chromSize[source_name] = source_size
            chain_id = None if len(fields) == 12 else fields[12] 

            if source_name not in maps:
                maps[source_name] = {} 
            if target_name not in maps[source_name]:
                maps[source_name][target_name] = {} 
            if source_start not in maps[source_name][target_name]: 
                maps[source_name][target_name][source_start] = {} 
            if target_start not in maps[source_name][target_name][source_start]: 
                maps[source_name][target_name][source_start][target_start] = [fields] 

            sfrom, tfrom = source_start, target_start 

        elif len(line.split()) == 3: 
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2]) 
            maps[source_name][target_name][source_start][target_start].append([int(fields[0]), int(fields[1]), int(fields[2])]) 

        elif len(line.split()) == 1: 
            maps[source_name][target_name][source_start][target_start].append([int(fields[0])]) 

    return maps 


ce_old_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta")) 
ce_new_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta")) 

chain = read_chain_file(sys.argv[4]) 

if verbose: 
    print(ce_old_dict.keys()) 
    print(ce_new_dict.keys()) 


if sys.argv[1] == "verify": 

    no_mismatches = True 
    for source in chain:
        #if (verbose): print(source) 
        if source not in ce_old_dict: 
            continue 
        for target in chain[source]: 
            #if (verbose): print(target) 
            if target not in ce_new_dict:
                continue 

            for source_start in chain[source][target]:
                sfrom = source_start 
                for target_start in chain[source][target][source_start]: 
                    fields = chain[source][target][source_start][target_start][0] 
                    
                    tfrom = target_start 
                    
                    if fields[9] == "-": 
                        tfrom = int(fields[8]) - int(fields[10]) - 1 

                    for el in chain[source][target][source_start][target_start][1:]: 
                        size = 0
                        sgap = 0
                        tgap = 0 
                        if len(el) == 3: 
                            size = el[0] 
                            sgap = el[1]
                            tgap = el[2] 
                        elif len(el) == 1: 
                            size = el[0] 

                        if len(ce_old_dict[source]) < sfrom+size+1:
                            break 

                        if fields[9] == "+": 
                            if len(ce_new_dict[target]) < tfrom+size+1:
                                break  
                        # TODO bounds checking for negative strand..     

                        
                        #print(ce_old_dict[source])
                        #print(ce_old_dict[source][sfrom:sfrom+size+1]) 
                        this_block_errors = 0 
                        dists = dict() 
                        last_idx_error = -1 
                        match = True

                        #print(fields)
                        #print("sfrom: " + str(sfrom)) 
                        #print("tfrom: " + str(tfrom)) 
                        #print("size:  " + str(size)) 
                        for idx in range(0, size): 
                            if (fields[9] == "+" and ce_old_dict[source][sfrom+idx].lower() != ce_new_dict[target][tfrom+idx].lower()) or (fields[9] == "-" and ce_old_dict[source][sfrom+idx].lower() != rev_comp[ce_new_dict[target][tfrom-idx].lower()]): 
                                #print("idx: " + str(idx) + " old: " + str(ce_old_dict[source][sfrom+idx].lower()) + " new: " + rev_comp[ce_new_dict[target][tfrom-idx].lower()])
                                match = False 
                                this_block_errors += 1  
                                if last_idx_error == -1:
                                    last_idx_error = idx 
                                else: 
                                    if idx - last_idx_error not in dists: 
                                        dists[idx - last_idx_error] = 1 
                                    else:
                                        dists[idx - last_idx_error] += 1 
                                        
                        #small_dist_count = 0 
                        #for el in dists.keys():
                        #    if el < 100:
                        #        small_dist_count += dists[el] 
                        #    else: 
                        #        break 

                        if (not match): 
                            #print("small_dist_count: " + str(small_dist_count)) 
                            print("mismatch: " + str(this_block_errors) + " errors in " + str(size) + " block" + "target: " + target + " header: " + " ".join(chain[source][target][source_start][target_start][0])) 
                            no_mismatches = False 
                        #else:
                            #print ("match") 

                        sfrom += size + sgap 
                        if (fields[9] == "+"): 
                            tfrom += (size + tgap) 
                        else: 
                            tfrom -= (size + tgap) 
    if no_mismatches: 
        print("PASS VERIFICATION") 


elif sys.argv[1] == "build": 
    new_chain = [] 
    if (verbose): print("") 
    if (verbose): print("=============\nFIRST PASS\n============") 
    for source in chain: 
        if source not in ce_old_dict: 
            print("source not in ce_old_dict") 
            exit(0) 

        for target in chain[source]: 
            if target not in ce_new_dict: 
                print("target not in ce_new_dict") 
                exit(0) 

            for source_start in chain[source][target]: 
                sfrom = source_start 
                for target_start in chain[source][target][source_start]: 
                    fields = chain[source][target][source_start][target_start][0] 
                    tfrom = target_start
                    
                    if fields[9] == "-": 
                        tfrom = int(fields[8]) - int(fields[10]) - 1 

                    if (verbose): print("\t".join(fields)) 
                    new_chain.append("\t".join(fields)) 
                    header_printed = False 
                    for chidx, el in enumerate(chain[source][target][source_start][target_start][1:]): 
                        #print(el) 
                        size = 0
                        sgap = 0
                        tgap = 0 
                        if len(el) == 3: 
                            size = el[0] 
                            sgap = el[1]
                            tgap = el[2] 
                        elif len(el) == 1: 
                            size = el[0] 

                        #print("sfrom+size: " + str(sfrom + size)) 
                        #print("tfrom+size: " + str(tfrom + size)) 
                        if len(ce_old_dict[source]) < sfrom+size:
                            print("lengths do not match")
                            exit(0) 

                        if fields[9] == "+": 
                            if len(ce_new_dict[target]) < tfrom+size:
                                print("lengths do not match")
                                exit(0) 
                        # TODO boundary check for "-" 

                        
                        #print(ce_old_dict[source])
                        #print(ce_old_dict[source][sfrom:sfrom+size+1]) 
                        if (fields[9] == "+" and not compare_seqs(ce_old_dict[source][sfrom:sfrom+size+1], ce_new_dict[target][tfrom:tfrom+size+1])) or (fields[9] == "-" and not compare_revcomp_seqs(ce_old_dict[source][sfrom:sfrom+size+1], ce_new_dict[target][tfrom - (size-1):tfrom+1])): 
                            #if ce_old_dict[source][sfrom:sfrom+size+1] != ce_new_dict[target][tfrom:tfrom+size+1]:
                            this_block_errors = 0 
                            last_idx_error = -1 
                            done = False 
                            for idx in range(0, size): 
                                # not same 
                                if (fields[9] == "+" and ce_old_dict[source][sfrom+idx].lower() != ce_new_dict[target][tfrom+idx].lower()) or (fields[9] == "-" and ce_old_dict[source][sfrom+idx].lower() != rev_comp[ce_new_dict[target][tfrom-idx].lower()]): 
                                    #if (verbose): print("broken idx: " + str(idx) + " " + str(ce_old_dict[source][sfrom+idx]) + " " + str(ce_new_dict[target][tfrom+idx])) 
                                    this_block_errors += 1  
                                    new_size = 0 
                                    if last_idx_error == -1:
                                        new_size = idx 
                                    else: 
                                        new_size = idx - last_idx_error - 1

                                    if len(el) == 3: 
                                        if (idx == size): 
                                            this_line = str(new_size) + " " + str(1 + sgap) + " " + str(1 + tgap) 
                                            if (verbose): print(this_line) 
                                            new_chain.append(this_line) 
                                            done = True 
                                        else: 
                                            this_line = str(new_size) + " 1 1" 
                                            if (verbose): print(this_line) 
                                            new_chain.append(this_line) 

                                    elif len(el) == 1: 
                                        if (idx == size): 
                                            if (verbose): print(str(new_size)) 
                                            new_chain.append(str(new_size)) 
                                            done = True 
                                        else: 
                                            this_line = str(new_size) + " 1 1" 
                                            if (verbose): print(this_line) 
                                            new_chain.append(this_line) 

                                    last_idx_error = idx 

                            if not done: 
                                if len(el) == 3: 
                                    this_line = str(size - last_idx_error - 1) + " " + str(sgap) + " " + str(tgap) 
                                    if (verbose): print(this_line) 
                                    new_chain.append(this_line) 
                                elif len(el) == 1: 
                                    this_line = str(size - last_idx_error - 1) 
                                    if (verbose): print(this_line) 
                                    new_chain.append(this_line) 
                                            
                            #print("mismatch: " + str(this_block_errors) + " errors in " + str(size) + " block") 
                        else:
                            #print ("match") 
                            if len(el) == 1: 
                                if (verbose): print(size) 
                                new_chain.append(str(size)) 
                            elif len(el) == 3: 
                                this_line = " ".join([str(x) for x in el])
                                if (verbose): print(this_line) 
                                new_chain.append(this_line) 


                        sfrom += size + sgap 
                        if (fields[9] == "+"): 
                            tfrom += size + tgap
                        else: 
                            tfrom -= (size + tgap) 



    new_chain2 = [] 
    last_el = "" 
    new_sgap = 0
    new_tgap = 0 
    last_size = -1 
    #######################################################
    # second pass analysis. remove all the tiny holes.. 
    ######################################################
    if (verbose): print("===============\nSECOND PASS\n=============") 
    for el in new_chain: 
        fields = el.split() 
        if "chain" in fields:  
            if (verbose): print(el) 
            new_chain2.append(el) 
            last_size = -1 
            new_sgap = 0
            new_tgap = 0 
        elif len(fields) == 3: 
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2]) 

        elif len(fields) == 1: 
            size = int(fields[0]) 


        if len(fields) == 3: 
            if size < min_region_size: 
                new_sgap += sgap + size 
                new_tgap += tgap + size 
                if (last_size == -1):  
                    last_size = 0 
            else: 
                if (last_size != -1): 
                    this_line = str(last_size) + "\t" + str(new_sgap) + "\t" + str(new_tgap)
                    if (verbose): print(this_line) 
                    new_chain2.append(this_line) 
                new_sgap = sgap  
                new_tgap = tgap 
                last_size = size 

        if len(fields) == 1: 
            if (last_size != -1): 
                this_line = str(last_size) + "\t" + str(new_sgap) + "\t" + str(new_tgap) 
                if (verbose): print(this_line) 
                new_chain2.append(this_line) 
            new_chain2.append(str(size)) 
            if (verbose): print(size) 
            #last_size = size 
                
        last_el = el 
    #if (verbose): print(str(last_size)) 
    #new_chain2.append(str(last_size)) 

    new_chain3 = [] 
    last_chain_header_idx = 0 

    ############################################
    # third pass analysis: update the headers 
    ############################################
    new_idx = 0
    this_strand = "+" 
    for idx,el in enumerate(new_chain2): 
        if (verbose): print(el) 
        fields = el.split() 
        if "chain" in fields: 
            last_chain_header_idx = new_idx 
            new_chain3.append(el) 
            new_idx += 1 
            this_strand = fields[9] 
        elif len(fields) == 3: 
            if fields[0] == "0": 
                # just change last header info 
                s_size = int(fields[1]) 
                t_size = int(fields[2]) 
                # add s_size to header start location locs 5 
                # add t_size to header start location loc 10 
                
                t_start = [int(new_chain3[last_chain_header_idx].split("\t")[10])]
                t_end = [int(new_chain3[last_chain_header_idx].split("\t")[11])]
                t_start = [str(t_start[0] + t_size)]

                new_chain3[last_chain_header_idx] = "\t".join(new_chain3[last_chain_header_idx].split("\t")[0:5] + [str(int(new_chain3[last_chain_header_idx].split("\t")[5]) + s_size)] + new_chain3[last_chain_header_idx].split("\t")[6:10] + [str(t_start[0])] + [str(t_end[0])] + new_chain3[last_chain_header_idx].split("\t")[12:]) 
            else:
                new_chain3.append(el) 
                new_idx += 1
                
        elif len(fields) == 1: 
            if fields[0] == "0": 
                # just change last header info 
                if ("chain" not in new_chain3[-1]): 
                    # must be a 3 field entry. 
                    s_size = int(new_chain3[-1].split()[1]) 
                    t_size = int(new_chain3[-1].split()[2]) 
                    
                    # remove last_size from header end locations  locs 6 and 11 
                    t_start = [int(new_chain3[last_chain_header_idx].split("\t")[10])]
                    t_end = [int(new_chain3[last_chain_header_idx].split("\t")[11])]
                    #t_end = [str(int(new_chain3[last_chain_header_idx].split("\t")[11]) - last_size)]
                    t_end = [str(t_end[0] - t_size)] 

                    new_chain3[last_chain_header_idx] = "\t".join(new_chain3[last_chain_header_idx].split("\t")[0:6] + [str(int(new_chain3[last_chain_header_idx].split("\t")[6]) - s_size)] + new_chain3[last_chain_header_idx].split("\t")[7:10] + [str(t_start[0])] + [str(t_end[0])] + new_chain3[last_chain_header_idx].split("\t")[12:]) 
                    new_chain3[-1] = new_chain3[-1].split()[0] 
                
            else:
                new_chain3.append(el) 
                new_idx += 1 

    if (verbose): print("===============\nTHIRD PASS\n=============") 
    for idx, el in enumerate(new_chain3): 
        if "chain" in el and idx > 0: 
            print() 
        print(el) 
