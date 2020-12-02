# python gaps_to_fasta.py <extracted gaps from references> <anticipated read size> <output fasta filename> <skip_length> 

import sys
import os.path 
 
if not len(sys.argv) == 5: 
    print("Not enough arguments.") 
    print("Usage: python gaps_to_fasta.py <extracted gaps from references> <anticipated read size> <output fasta filename> <skip_length>") 
    exit() 

in_fn = sys.argv[1] 
read_size = int(sys.argv[2]) 
out_fn = sys.argv[3] 
skip_len = int(sys.argv[4])
if skip_len < 1: 
    print("not a valid skip length" )
    exit() 

if not os.path.isfile(in_fn): 
    print("not a file: " + in_fn )
    exit() 

if os.path.isfile(out_fn): 
    prompt = "overwrite " + sys.argv[3] + "? [Enter] for yes [any other key] for no" 
    text = input(prompt) 
    if text != "":
        exit() 

outfile = open(out_fn, "w") 

seq_name = "" 
seq_num = 0 
seq = "" 
for line in open(in_fn): 
    if ">" in line: 
        if len(seq) > 0: # if we have a sequence, add the seqs to our fasta file. 
            for i in range(read_size, len(seq), skip_len): 
                outfile.write(">" + seq_name.strip() + "_" + str(seq_num) + "\n") # sequence name 
                outfile.write(seq[seq_num:i] + "\n") # actual sequence 
#                print(seq[seq_num:i] + "\n") 
                seq_num = seq_num + skip_len 
        seq_name = line.split(">")[1] 
        seq_num = 0 
        seq = "" 
    else: 
        seq = seq + line.strip()
        

if len(seq) > 0: # if we have a sequence, add the seqs to our fasta file. 
    for i in range(read_size, len(seq), skip_len): 
        outfile.write(">" + seq_name.strip() + "_" + str(seq_num) + "\n") # sequence name 
        outfile.write(seq[seq_num:i] + "\n") # actual sequence 
#        print(seq[seq_num:i] + "\n") 
        seq_num = seq_num + skip_len 
seq_name = line.split(">")[1] 
seq_num = 0 
seq = "" 
