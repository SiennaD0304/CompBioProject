##### INPUT: list of motifs, sequence alignment results
##### OUTPUT: correctly formatted file to run MEME suite

import pandas as pd

#### FILTER MOTIFS AND HITS

def make_filtered_dict(in_file):

    # make df of input file, extract info needed
    with open(in_file, 'r') as f:
        df = pd.read_csv(f, sep='\t')
        motifs=set(df['Motif'].tolist()) # only want unique motifs 
        aligned_motifs = df['Aligned_Motif'].tolist()
        aligned_hits = df['Aligned_Hit'].tolist()

        #dictionary with aligned motifs as keys and aligned hits as values
        motif_hit_dict={}
        for i in range(len(aligned_motifs)):
            motif_hit_dict[aligned_motifs[i]]=aligned_hits[i]

        # filtering by length
        filtered_dict={}
        for key, value in motif_hit_dict.items():
            if len(key)>=4:
                filtered_dict[key]=value

    return motifs, filtered_dict

#### DESIGN PROBABILITY MATRIX

def make_prob_matrix(motif, filtered_dict):
    alphabet = "ACDEFGHIKLMNPQRSTVWY" # amino acids
    
    # initialize count matrix to zeros for specific motif length
    count_matrix = pd.DataFrame(0, index=range(len(motif)), columns=list(alphabet))

    # initialize with the motif itself
    for i, aa in enumerate(motif):
        if aa in count_matrix.columns: # match positions in motif to amino acids in alphabet
            count_matrix.loc[i, aa] += 1

    # update matrix based on substitutions in aligned hits
    for fragment, hit in filtered_dict.items():
        start = motif.find(fragment) # find where fragment is in the entire motif
        if start != -1: # check to make sure there's a match 
            for i, char in enumerate(hit):
                position = start + i 
                if position < len(motif) and char in alphabet:
                    count_matrix.loc[position, char] += 1 # update positions with an aligned hit
        
    # normalize: divide each position by the row total 
    row_sums = count_matrix.sum(axis=1)
    prob_matrix = count_matrix.divide(row_sums, axis=0).round(3)
        
    return prob_matrix

#### WRITE MEME INPUT FILE

def write_motif_file(motifs, filtered_dict, out_file):
    alphabet = "ACDEFGHIKLMNPQRSTVWY" # amino acids

    with open(out_file, "w") as f:
        # write header
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        # information for each specific motif
        for i, motif in enumerate(motifs):
            f.write(f"MOTIF {motif} motif_{i+1}\n")
            
            # call probability matrix function for each dataset
            probs = make_prob_matrix(motif, filtered_dict)
            
            # write motif metadata and probability matrix to file
            f.write(f"letter-probability matrix: alength= {len(alphabet)} w= {len(motif)} nsites= 1\n")
            f.write(f"{probs.to_string(header=False, index=False)}\n\n")

    return f

# process input file
motifs_set, filtered_data = make_filtered_dict("msa_results.tsv")
# write output file
write_motif_file(motifs_set, filtered_data, "motifs_input.meme")