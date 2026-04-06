from Bio import SeqIO
###############################
# Used the Module Pymemesuite which allows you to access meme suit's tool FIMO:Find Individual Motif Occurrences
    #FIMO: scans a set of sequences for individual matches to motifs you provide 
################################
#LINK TO SITE WITH INSTALLION INSTRUCTIONS FOR MODULE
#https://pypi.org/project/pymemesuite/#:~:text=%F0%9F%94%A7%20Installing,qvalue%20)
# from command line installed the module pymemesuite with the command:
#pip install pymemesuite 

import requests 
import Bio.SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

####################
# GET MUSCLE PROTEIN SEQUENCES
#####################
proteins = {
    "ABLIM1": "O14639",
    "MYBPC": "Q14896",
    "MYL2": "P10916"
}

def fetch_fasta(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    return requests.get(url).text

sequences = {}
fasta_list = []

for name, uid in proteins.items():
    fasta = fetch_fasta(uid)
    fasta_list.append(fasta)
    record = SeqIO.read(StringIO(fasta), "fasta")
    sequences[name] = str(record.seq)
    #seq_list.append(sequences[name])

fasta_files = "/n".join(fasta_list)
fasta_input = fasta_files.replace("/n", "")

with open ("muscle_seq.fasta", "w") as f: 
    f.write(fasta_input)

# given fasta files and list of motifs
# generate correct format for inputs

# first, need to update headers of muscle protein fasta file

# function to format input for motifs
# edit to include background frequency lines
def write_motif_file(motifs, out_file):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    
    with open(out_file, "w") as f:
        # write header
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, residue in enumerate(motifs):
            # motif identifier and alternate name
            f.write(f"MOTIF {residue} motif_{i+1}\n")

            # letter-probability matrix header
            f.write(f"letter-probability matrix: alength= {len(alphabet)} w= {len(residue)} nsites= 1\n")
            
            # design probability matrix
            # shows the probabilities of each amino acid at each position in motif
            # this is a strict matrix - edit to account for subsitutions
            for i in residue:
                row = []
                for aa in alphabet:
                    if i == aa:
                        row.append("1.0000")
                    else:
                        row.append("0.0000")
                f.write("  " + "  ".join(row) + "\n")

# test 
motifs = ["GTVLSDIIQKYFF", "LSDLIQ", "LSSLIQ", "LSDIIQ", "GLALSDLIQKYFF"]
write_motif_file(motifs, "motifs_input.meme")


#### START OF PYMEMESUITE TOOL #### 
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
    for record in Bio.SeqIO.parse("muscle_seq.fasta", "fasta")
]

# dictionary to store ids and sequences
seq_id_dict = {}
for record in Bio.SeqIO.parse("muscle_seq.fasta", "fasta"):
    seqs = str(record.seq)
    name = record.id
    seq_id_dict.update({name:seqs}) 

fimo = FIMO(both_strands=False)

sequence_pts = {}
pdb_names_ids = []
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


                temp_list= m.start, m.stop
                sequence_pts.update({match_name:temp_list})
                sequence_keys = list(sequence_pts.keys())

                for motifs in sequence_keys: 
                #for ABLIM1 there is no pdb number and therefore must be downloaded before from AlphaFold 
                #To download go to: https://alphafold.ebi.ac.uk/entry/AF-O14639-4-F1
                #click download and select pdb 
                #save to preferred file 
                #open PyMOL and type command: 
                #load (your/path/filename.pdb)

                    if motifs == "sp|O14639|ABLM1_HUMAN": 
                        #pdb_names_ids.append("https://rest.uniprot.org/uniprotkb/O14639.fasta")
                        pymol_script_1 = ["load /Users/siennad0304/Documents/Senior Year/Comp Biology /AF-O14639-F1-model_v6.pdb",
                                          "as cartoon", 
                                          (f"select i. {m.start}-{m.stop}"), 
                                          "set_name sele, motif",
                                          "hide everything",
                                          "show sticks, motif",
                                          "fetch 6zu0", #this this the Nav1.5 sodium channel which we will align to  
                                          "align 6zu0, motif",
                                          "orient motif"

                                            ]
                        with open(f"{motifs}.pml", "w") as l:
                            l.write("\n".join(pymol_script_1))

                    if motifs == "sp|Q14896|MYPC3_HUMAN":
                        #pdb_names_ids.append("1PD6")
                        pymol_script_2 = ["fetch 1PD6",
                                          "as cartoon", 
                                          (f"select i. {m.start}-{m.stop}"), 
                                          "set_name sele, motif",
                                          "hide everything",
                                          "show sticks, motif",
                                          "fetch 6zu0", #this this the Nav1.5 sodium channel which we will align to  
                                          "align 6zu0, motif",
                                          "orient motif"
                                            ]
                        with open(f"{motifs}.pml", "w") as x:
                            x.write("\n".join(pymol_script_2))

                    if motifs == "sp|P10916|MLRV_HUMAN":
                        #pdb_names_ids.append("5TBY")
                        pymol_script_3 = ["fetch 5TBY"
                                          "as cartoon", 
                                          (f"select i. {m.start}-{m.stop}"), 
                                          "set_name sele, motif",
                                          "hide everything",
                                          "show sticks, motif",
                                          "fetch 6zu0", #this this the Nav1.5 sodium channel which we will align to  
                                          "align 6zu0, motif",
                                          "orient motif"
                                            ]
                        
                        with open(f"{motifs}.pml", "w") as z:
                            z.write("\n".join(pymol_script_3))
