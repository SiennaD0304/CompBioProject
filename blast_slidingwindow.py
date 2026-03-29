###################
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

####################
# GET SEQUENCES
#####################
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
motifs = ["GLALSDLIQKYFF", "LSDLIQ", "LSSLIQ"]


############################
# SLIDING WINDOW ATTEMPT
############################
#score match function
def score_match(seq1, seq2):
    score = 0
    for a, b in zip(seq1, seq2):
        if a == b:
            score += 2
        elif a in "ST" and b in "ST":
            score += 1
        else:
            score -= 1
    return score

results = []

#scan through proteins
for protein in ["ABLIM1", "MYBPC", "MYL2"]:
    seq = sequences[protein]
    for motif in motifs:
        m_len = len(motif)
        for i in range(len(seq) - m_len + 1):
            window = seq[i:i+m_len]
            score = score_match(window, motif)

            if score > 5: #threshold
                results.append({
                    "Protein": protein,
                    "Start": i,
                    "End": i + m_len,
                    "Motif": motif,
                    "Match": window,
                    "Score": score
                    })

#put results in pandas               
df = pd.DataFrame(results)
df = df.sort_values(by="Score", ascending=False)

print(df.head(20))

df.to_csv("motif_scan_results.tsv", sep = "\t")


#############################
#BLAST ATTEMPT
#############################

SeqIO.write(
    SeqRecord(Seq(sequences["Nav1.5"]), id="Nav1.5"),
    "Nav1.5.fasta",
    "fasta"
)

muscle_records = [
    SeqRecord(Seq(sequences["ABLIM1"]), id="ABLIM1"),
    SeqRecord(Seq(sequences["MYBPC"]), id="MYBPC"),
    SeqRecord(Seq(sequences["MYL2"]), id="MYL2"),
]

SeqIO.write(muscle_records, "muscle_proteins.fasta", "fasta")

motif_positions = []
nav_seq = sequences["Nav1.5"]

for motif in motifs:
    start_idx = nav_seq.find(motif)
    if start_idx != -1:
        motif_positions.append({
            "Motif": motif,
            "Start": start_idx + 1,  # +1 for 1-based indexing like BLAST
            "End": start_idx + len(motif)
        })

#create blast db
os.system("makeblastdb -in muscle_proteins.fasta -dbtype prot -out muscle_db")
#actual blast command
os.system("blastp -query Nav1.5.fasta -db muscle_db -out nav_vs_muscle.tsv -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'")

#put blast results into pandas for potential later use
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

#attempt to add specific motifs and matches to blast results, but blast doesn't really do motifs?
def assign_motif(row):
    for mp in motif_positions:
        #checking if motif overlaps query alignment
        if (row['Query Start'] <= mp['End']) and (row['Query End'] >= mp['Start']):
            return mp['Motif']
    return "No motif"

blast_df["Motif"] = blast_df.apply(assign_motif, axis=1)

def matched_seq(row):
    subj_seq = sequences[row["Subject"]]
    start = row["Subject Start"] - 1  # convert to 0-based
    end = row["Subject End"]
    return subj_seq[start:end]

blast_df["Motif Match"] = blast_df.apply(matched_seq, axis=1)

print("\nBLAST Results:")
print(blast_df.head())

blast_df.to_csv("nav_vs_muscle.tsv", sep="\t", index=False)

#############################
#MULTIPLE SEQUENCE ALIGNMENT
#aligning each BLAST hit region against all three Nav1.5 motifs
#using biopython pairwise aligner with blosum62
#output one row per pair (blast hit + nav1.5 motif), showing alignment stats
#############################

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "global"

blosum62 = substitution_matrices.load("BLOSUM62")
msa_rows = []

for _, row in blast_df.iterrows():
    subject_name = row["Subject"]
    blast_hit_seq = row["Motif Match"]

    for motif in motifs:
        alignment = next(iter(aligner.align(motif, blast_hit_seq)))

        aligned_motif = alignment[0]
        aligned_hit = alignment[1]

        #calculate percent identity- matching residues / aligned (non-gap) positions
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
            "Alignment_Score": round(alignment.score, 1),
        })

msa_df = pd.DataFrame(msa_rows)

print("\nMSA Results:")
print(msa_df.to_string(index=False))

msa_df.to_csv("msa_results.tsv", sep="\t", index=False)
print("\nMSA results written to: msa_results.tsv")


###### TO RUN: 
###in terminal
# cd CompBioProject
# python blast_slidingwindow.py