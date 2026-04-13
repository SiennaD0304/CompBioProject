#the file name will need to be changed to your path before running
#for ABLIM1 there is no pdb number and therefore must be downloaded before from AlphaFold 
#to download file:https://alphafold.ebi.ac.uk/entry/AF-O14639-F1
#click download and select pdb 
#save to preferred file 
load /Users/siennad0304/Documents/Senior Year/Comp Biology /AF-O14639-F1-model_v6.pdb
as cartoon
select i. 126-138
set_name sele, motif
hide everything
show sticks, motif
fetch 6zu0
align 6zu0, motif
orient motif
