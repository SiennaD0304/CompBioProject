# Pymemesuite module to access meme suite's FIMO tool
# FIMO scans set of sequences for matches to provided motifs

# MODULE DOCUMENTATION
#https://pypi.org/project/pymemesuite/#:~:text=%F0%9F%94%A7%20Installing,qvalue%20)

# from virtual environment on command line installed the module pymemesuite with the command:
# pip install pymemesuite 

import Bio.SeqIO
from Bio import SeqIO
from Bio.Seq import Seq

############################
#START OF MEME SUITE MODULE:FIMO MOTIF SEARCH 
############################

from pymemesuite.common import MotifFile
from pymemesuite.common import Sequence 
from pymemesuite.fimo import FIMO 

# biopython to extract amino acid sequences of cardiac muscle proteins
sequences = [
    Sequence(str(record.seq), name=record.id.encode())
    for record in Bio.SeqIO.parse("muscle_seq.fna", "fasta")
]

# dictionary to store ids and sequences
seq_id_dict = {}
for record in Bio.SeqIO.parse("muscle_seq.fna", "fasta"):
    seqs = str(record.seq)
    name = record.id
    seq_id_dict.update({name:seqs}) 

fimo = FIMO(both_strands=False)

# WRITE RESULTS
with open("meme_suite_output.txt", "w") as out:

    # open input file to iterate through all motifs
    # MotifFile function part of pymemesuite.common
    # this loads meme format file without manual parsing
    with MotifFile("motifs_input.meme") as motif_file:

        for motif in motif_file: # iterate through motifs
            out.write(f"{motif.name.decode()} {motif.consensus} \n")
            header = "Accession\tstart\tstop\tstrand\tscore\tpvalue\tmatch_seq\n"
            out.write(header)

            pattern = fimo.score_motif(motif, sequences, motif_file.background) # run FIMO

            for m in pattern.matched_elements:
                match_name = m.source.accession.decode() # prot with motif match
                full_match_seq = seq_id_dict[match_name] # sequence with motif match

                sub_seq = full_match_seq[m.start-1 : m.stop] # this is the actual motif match! need to splice full sequence

                line = f"{match_name}\t{m.start}\t{m.stop}\t{m.strand}\t{m.score}\t{m.pvalue}\t{sub_seq}\n"
                out.write(line)
            
            out.write("\n")