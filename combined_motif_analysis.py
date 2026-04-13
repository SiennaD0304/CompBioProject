##################
# IMPORT PACKAGES
##################
import requests
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio.Align import substitution_matrices
from io import StringIO
import pandas as pd
import os
from pymemesuite.common import MotifFile
from pymemesuite.common import Sequence
from pymemesuite.fimo import FIMO

##################
# GET PROTEINS
##################
proteins = {
    "Nav1.5": "P15389",
    "ABLIM1": "O14639",
    "MYBPC": "Q14896",
    "MYL2": "P10916"
}

def fetch_fasta(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    return requests.get(url).text

sequences = {}

for name, uid in proteins.items():
    fasta = fetch_fasta(uid)
    record = SeqIO.read(StringIO(fasta), "fasta")
    sequences[name] = str(record.seq)

#define motifs
motifs = ["GLALSDLIQKYFF", "LSDLIQ", "LSSLIQ", "GTVLSDIIQKYFF", "LSDIIQ", "FKVGHGLAC", "VNILAKI"]

###################
# BLAST
################
SeqIO.write(
    SeqRecord(Seq(sequences["Nav1.5"]), id = "Nav1.5"),
    "Nav1.5.fasta",
    "fasta"
)

muscle_records = [
    SeqRecord(Seq(sequences["ABLIM1"]), id = "ABLIM1"),
    SeqRecord(Seq(sequences["MYBPC"]), id = "MYBPC"),
    SeqRecord(Seq(sequences["MYL2"]), id = "MYL2")
]

SeqIO.write(muscle_records, "muscle_proteins.fasta", "fasta")

motif_positions = []
nav_seq = sequences["Nav1.5"]

for motif in motifs:
    start_idx = nav_seq.find(motif)
    if start_idx != -1:
        motif_positions.append({
            "Motif": motif,
            "Start": start_idx + 1,
            "End": start_idx + len(motif)
        })

os.system("makeblastdb -in muscle_proteins.fasta -dbtype prot -out muscle_db")

os.system("blastp -query Nav1.5.fasta -db muscle_db -out nav_vs_muscle.tsv -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'")

blast_df = pd.read_csv(
    "nav_vs_muscle.tsv",
    sep = "\t",
    names = [
        "Query",
        "Subject",
        "Percent Identity",
        "Alignment Length",
        "Query Start",
        "Query End",
        "Subject Start",
        "Subject End",
        "E-value",
        "Bit Score"
    ]
)

#adding specific motifs/matches to blast results
def assign_motif(row):
    for mp in motif_positions:
        #check if motif overlaps query alignment
        if (row['Query Start'] <= mp['End']) and (row['Query End'] >= mp['Start']):
            return mp['Motif']
    return "No motif"

blast_df["Motif"] = blast_df.apply(assign_motif, axis=1)

def matched_seq(row):
    subj_seq = sequences[row["Subject"]]
    start = row["Subject Start"] - 1 #cinvert to 0 based
    end = row["Subject End"]
    return subj_seq[start:end]

blast_df["Motif Match"] = blast_df.apply(matched_seq, axis=1)
blast_df.to_csv("nav_vs_muscle_results.tsv", sep="\t", index=False)

#####################
# ALIGNMENT
#####################
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "local"

blosum62 = substitution_matrices.load("BLOSUM62")
msa_rows = []

for _, row in blast_df.iterrows():
    subject_name = row["Subject"]
    blast_hit_seq = row["Motif Match"]

    for motif in motifs:
        alignment = next(iter(aligner.align(motif, blast_hit_seq)))
        
        aligned_motif = alignment[0]
        aligned_hit = alignment[1]

        #calculate percent identity
        comparable = sum(1 for a, b in zip(aligned_motif, aligned_hit) if a != "-" and b != "-")
        identities = sum(a == b for a, b in zip(aligned_motif, aligned_hit) if a != "-" and b != "-")
        similarities = sum(
            1 for a, b in zip(aligned_motif, aligned_hit)
            if a != "-" and b != "-" and blosum62.get((a, b), -99) > 0
        )
        pct_identity = round(100.0 * identities / comparable, 1) if comparable else 0.0
        pct_similarity = round(100.0 * similarities / comparable, 1) if comparable else 0.0

        msa_rows.append({
            "Subject": subject_name,
            "BLAST_Hit_Region": blast_hit_seq,
            "Motif": motif,
            "Aligned_Motif": aligned_motif,
            "Aligned_Hit": aligned_hit,
            "Identity(%)": pct_identity,
            "Similarity(%)": pct_similarity,
            "Alignment_Score": round(alignment.score, 1)
        })

msa_df = pd.DataFrame(msa_rows)
msa_df.to_csv("msa_results.tsv", sep="\t", index=False)

########################
#CREATING MEME INPUT
######################

# file 1
def write_motif_file(motif, out_file):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    with open(out_file, "w") as f:
        f.write("MEME version 5\n")
        f.write(f"ALPHABET = {alphabet}\n\n")
        for i, residue in enumerate(motifs):
            f.write(f"MOTIF {residue} motif_{i+1}\n")
            f.write(f"letter-probability matrix: alength = {len(alphabet)} w = {len(residue)} nsites = 1\n")

            for i in residue: # strict matrix
                row = []
                for aa in alphabet:
                    if i == aa:
                        row.append("1.0000")
                    else:
                        row.append("0.000")
                f.write("  " + "  ".join(row) + "\n")

write_motif_file(motifs, "motifs_strict_input.meme")

# file 2
def make_filtered_dict(in_file):

    with open(in_file, 'r') as f:
        df = pd.read_csv(f, sep='\t')
        motifs=set(df['Motif'].tolist()) # only want unique motifs 
        aligned_motifs = df['Aligned_Motif'].tolist()
        aligned_hits = df['Aligned_Hit'].tolist()

        motif_hit_dict={}
        for i in range(len(aligned_motifs)):
            motif_hit_dict[aligned_motifs[i]]=aligned_hits[i]

        # filtering by length
        filtered_dict={}
        for key, value in motif_hit_dict.items():
            if len(key)>=4:
                filtered_dict[key]=value

    return motifs, filtered_dict

def make_prob_matrix(motif, filtered_dict):
    alphabet = "ACDEFGHIKLMNPQRSTVWY" 
    
    # initialize count matrix to zeros for specific motif length
    count_matrix = pd.DataFrame(0, index=range(len(motif)), columns=list(alphabet))

    # initialize with the motif itself
    for i, aa in enumerate(motif):
        if aa in count_matrix.columns: 
            count_matrix.loc[i, aa] += 1

    # update matrix based on substitutions in aligned hits
    for fragment, hit in filtered_dict.items():
        start = motif.find(fragment) 
        if start != -1: # check to make sure there's a match 
            for i, char in enumerate(hit):
                position = start + i 
                if position < len(motif) and char in alphabet:
                    count_matrix.loc[position, char] += 1 # update positions with an aligned hit

    row_sums = count_matrix.sum(axis=1)
    prob_matrix = count_matrix.divide(row_sums, axis=0).round(3) # normalize
        
    return prob_matrix

def write_flexible_motif_file(motifs, filtered_dict, out_file):
    alphabet = "ACDEFGHIKLMNPQRSTVWY" # amino acids

    with open(out_file, "w") as f:
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, motif in enumerate(motifs):
            f.write(f"MOTIF {motif} motif_{i+8}\n")
            
            # call probability matrix function for each motif
            probs = make_prob_matrix(motif, filtered_dict)
            
            f.write(f"letter-probability matrix: alength= {len(alphabet)} w= {len(motif)} nsites= 1\n")
            f.write(f"{probs.to_string(header=False, index=False)}\n\n")

    return f

motifs_set, filtered_data = make_filtered_dict("msa_results.tsv")
write_flexible_motif_file(motifs_set, filtered_data, "motifs_flexible_input.meme")


###############
#FIMO
###############

import shutil
shutil.copy("muscle_proteins.fasta", "muscle_seq.fna")

sequences_meme = [
    Sequence(str(record.seq), name=record.id.encode())
    for record in SeqIO.parse("muscle_seq.fna", "fasta")
]

seq_id_dict = {}
for record in SeqIO.parse("muscle_seq.fna", "fasta"):
    seqs = str(record.seq)
    name = record.id
    seq_id_dict.update({name:seqs})

fimo = FIMO(both_strands=False)

meme_results = []
with open("meme_suite_output.txt", "w") as out:
    with MotifFile("motifs_strict_input.meme") as motif_file:
        for motif in motif_file:
            out.write(f"{motif.name.decode()} {motif.consensus} \n")
            header = "Accession\tstart\tstop\tstrand\tscore\tpvalue\tmatch_seq\n"
            out.write(header)
            pattern = fimo.score_motif(motif, sequences_meme, motif_file.background)
            for m in pattern.matched_elements:
                match_name = m.source.accession.decode()
                full_match_seq = seq_id_dict[match_name]
                sub_seq = full_match_seq[m.start-1 : m.stop]
                line = f"{match_name}\t{m.start}\t{m.stop}\t{m.strand}\t{m.score}\t{m.pvalue}\t{sub_seq}\n"
                out.write(line)

                meme_results.append({
                    "Method": "MEME",
                    "Protein": match_name,
                    "Motif": motif.name.decode(),
                    "Start": m.start,
                    "End": m.stop,
                    "Score": m.score,
                    "P_value": m.pvalue,
                    "Match_Sequence": sub_seq
                })
            
            out.write("\n")

    with MotifFile("motifs_flexible_input.meme") as motif_file2:
        for motif in motif_file2:
            out.write(f"{motif.name.decode()} {motif.consensus} \n")
            header = "Accession\tstart\tstop\tstrand\tscore\tpvalue\tmatch_seq\n"
            out.write(header)
            pattern = fimo.score_motif(motif, sequences_meme, motif_file2.background)
            for m in pattern.matched_elements:
                match_name = m.source.accession.decode()
                full_match_seq = seq_id_dict[match_name]
                sub_seq = full_match_seq[m.start-1 : m.stop]
                line = f"{match_name}\t{m.start}\t{m.stop}\t{m.strand}\t{m.score}\t{m.pvalue}\t{sub_seq}\n"
                out.write(line)

                meme_results.append({
                    "Method": "MEME",
                    "Protein": match_name,
                    "Motif": motif.name.decode(),
                    "Start": m.start,
                    "End": m.stop,
                    "Score": m.score,
                    "P_value": m.pvalue,
                    "Match_Sequence": sub_seq
                })
            
            out.write("\n")

print(f"Found {len(meme_results)} MEME matches")

#############
# COMBINED RESULTS ANALYSIS
#############

# Convert BLAST+MSA results to comparable format
blast_msa_results = []
for _, row in msa_df.iterrows():
    blast_msa_results.append({
        "Method": "BLAST+MSA",
        "Protein": row["Subject"],
        "Motif": row["Motif"],
        "Start": None,  # BLAST doesn't provide exact motif positions
        "End": None,
        "Score": row["Alignment_Score"],
        "Identity_Pct": row["Identity(%)"],
        "Similarity_Pct": row["Similarity(%)"],
        "Match_Sequence": row["Aligned_Hit"].replace("-", "")
    })

# Create comprehensive comparison DataFrame
all_results = []

# Add MEME results
for result in meme_results:
    all_results.append({
        "Method": result["Method"],
        "Protein": result["Protein"],
        "Motif": result["Motif"],
        "Start": result["Start"],
        "End": result["End"],
        "MEME_Score": result["Score"],
        "MEME_Pvalue": result["P_value"],
        "MSA_Score": None,
        "MSA_Identity": None,
        "MSA_Similarity": None,
        "Sliding_Score": None,
        "Match_Sequence": result["Match_Sequence"]
    })

# Add BLAST+MSA results
for result in blast_msa_results:
    all_results.append({
        "Method": result["Method"],
        "Protein": result["Protein"],
        "Motif": result["Motif"],
        "Start": result["Start"],
        "End": result["End"],
        "MEME_Score": None,
        "MEME_Pvalue": None,
        "MSA_Score": result["Score"],
        "MSA_Identity": result["Identity_Pct"],
        "MSA_Similarity": result["Similarity_Pct"],
        "Sliding_Score": None,
        "Match_Sequence": result["Match_Sequence"]
    })

combined_df = pd.DataFrame(all_results)

# Find best overall matches by looking for protein-motif combinations found by multiple methods
def identify_best_matches(df):
    best_matches = []
    
    # Group by protein and motif to find overlapping results
    grouped = df.groupby(["Protein", "Motif"])
    
    for (protein, motif), group in grouped:
        methods_used = group["Method"].tolist()
        
        # Calculate composite score based on available methods
        composite_score = 0
        score_components = []
        
        # MEME component (normalize p-value to score, lower p-value = higher score)
        meme_rows = group[group["Method"] == "MEME"]
        if len(meme_rows) > 0:
            meme_pval = meme_rows.iloc[0]["MEME_Pvalue"]
            if meme_pval is not None:
                meme_score = -np.log10(meme_pval) if meme_pval > 0 else 10
                composite_score += meme_score * 0.4
                score_components.append(f"MEME: {meme_score:.2f}")
        
        # MSA component
        msa_rows = group[group["Method"] == "BLAST+MSA"]
        if len(msa_rows) > 0:
            msa_score = msa_rows.iloc[0]["MSA_Score"]
            msa_identity = msa_rows.iloc[0]["MSA_Identity"]
            if msa_score is not None:
                composite_score += (msa_score/10) * 0.3  # Normalize MSA score
                score_components.append(f"MSA: {msa_score}")
            if msa_identity is not None:
                composite_score += (msa_identity/100) * 0.2  # Add identity component
                score_components.append(f"ID: {msa_identity}%")
        
        best_matches.append({
            "Protein": protein,
            "Motif": motif,
            "Methods_Used": ", ".join(methods_used),
            "Num_Methods": len(methods_used),
            "Composite_Score": round(composite_score, 3),
            "Score_Components": "; ".join(score_components),
            "Sequences": "; ".join([seq for seq in group["Match_Sequence"].unique() if pd.notna(seq)])
        })
    
    return pd.DataFrame(best_matches).sort_values("Composite_Score", ascending=False)

import numpy as np
best_matches_df = identify_best_matches(combined_df)

print(f"\nFound {len(best_matches_df)} unique protein-motif combinations")
print("\nTOP 10 BEST OVERALL MATCHES:")
print("="*80)
print(best_matches_df.head(10).to_string(index=False))

# Save all results
combined_df.to_csv("all_methods_combined.tsv", sep="\t", index=False)
best_matches_df.to_csv("best_overall_matches.tsv", sep="\t", index=False)

print(f"\n=== ANALYSIS COMPLETE ===")
print(f"Results saved to:")
print(f"- nav_vs_muscle_results.tsv (BLAST results)")
print(f"- msa_results.tsv (multiple sequence alignment)")
print(f"- meme_suite_output.txt (MEME results)")
print(f"- all_methods_combined.tsv (all results combined)")
print(f"- best_overall_matches.tsv (best matches ranked by composite score)")

##### TO RUN: 
###in terminal
# cd CompBioProject
# python combined_motif_analysis.py
