from tqdm import tqdm
import pandas as pd
from sequana import *
import numpy as np
import json

def entropy(sequence):
    pi = [sequence.count(l) / float(len(sequence)) for l in "ACGT"]
    pi = [x for x in pi if x != 0]
    return -sum(pi * np.log(pi))


def convert_all_numpy(obj):
    if isinstance(obj, dict):
        return {k: convert_all_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_all_numpy(v) for v in obj]
    elif isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64)):
        return float(obj)
    else:
        return obj

def save_json(filename, obj):
    import json
    obj = convert_all_numpy(obj)
    with open(filename, "w") as fout:
        json.dump(obj, fout)


def get_INDELs(bamfile, fasta_file, minimal_length_indels=100):
    
    data = []

    b1 = BAM(bamfile)
    N = len(BAM(bamfile).get_df().query("flag in [0,16]"))
    print(f"Number of mapped reads (flag 0,16) : {N}")
    b1 = BAM(bamfile)
    
    genome = FastA(fasta_file).sequences[0]

    N = 0
   
    for a in b1:
        if a.flag in [0,16]: #, 2048, 2064]:
            N+=1
            position_query = 0
            sequences = []
            position_ref = a.reference_start
            deletions = []
            for typ, length in Cigar(a.cigarstring).as_tuple():
                if typ=='S':
                    position_query += length
                    position_ref +=0
                elif typ == 'M':
                    sequences.append(('M', position_query, length, position_ref,a.query_sequence[position_query:position_query+length]))
                    position_query += length
                    position_ref += length
                elif typ == 'D':
                    sequences.append(('D', position_query, length, position_ref,genome[position_ref:position_ref+length]))
                    position_query += 0
                    # get the sequence deleted from the genome
                    position_ref += length
                elif typ == 'I':
                    sequences.append(('I', position_query, length, position_ref, a.query_sequence[position_query:position_query+length]))
                    position_query += length
                    position_ref += 0
                
            for typ, qpos, L, rpos, sequence in sequences:
                if typ=='I' and len(sequence)>minimal_length_indels:
                    data.append([typ,qpos,L,rpos,sequence, a.qname, a.rname])
                if typ=='D' and len(sequence)>minimal_length_indels:
                    data.append([typ,qpos,L,rpos,sequence, a.qname, a.rname])
    df = pd.DataFrame(data)
    if len(data):
        df.columns = ['typ', 'qpos', 'L', 'rpos', 'sequence', 'qname', 'rname']
        # total number of reads that map, N number of INDELs
        
        print("df length", len(df), f"N={N}", "%indels=", len(df) / N * 100)
        print("% of insertions", sum(df['typ'] == 'I') / N * 100, " Ninsertions", sum(df['typ'] == 'I')) 
        print("% of deletions", sum(df['typ'] == 'D') / N *100, " Ndeletions", sum(df['typ'] == 'D'))

        df['entropy'] = [entropy(x) for x in df.sequence.values]

        stats = {}
        stats['N'] = N
        stats['percent_insertions'] =  round(sum(df['typ'] == 'I') / N * 100,2)
        stats['percent_deletions'] =  round(sum(df['typ'] == 'D') / N * 100,2)
        stats['insertions'] =  sum(df['typ'] == 'I') 
        stats['deletions'] =  sum(df['typ'] == 'D')

        return df, stats
    else:
        stats = {
                'N': N, 
                'insertions': 0, 
                'deletions':0, 
                'percent_insertions': 0, 
                'percent_deletions': 0
            }
        return pd.DataFrame({'x':[], 'y':[]}), stats 


