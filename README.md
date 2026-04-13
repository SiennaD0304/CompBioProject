The entire analysis pipeline for this project is written in the combined_motif_analysis.py, which inputs a list of proteins and UniProt IDs to run BLAST, multiple sequence alignment, and MEME Suite FIMO to generate candidate motif matches. These candidate motifs are then used to generate PyMOL scripts to visualize each motif for putative interaction with the scorpion toxin LqhIII. 

This code requires a virtual environment. In VS Code, set up a virtual environment using command+shift+P to bring up the command palette. Then select "Python: Create Environment" and choose "Venv". VS Code will automatically create a .venv folder in your workspace and set it as the active environment.

Documentation on how to set up a virtual environment in VS Code can also be found here: https://code.visualstudio.com/docs/python/environments. 
Once the virtual environment is set up in your project folder, make sure that it is activated using the command:
```
source .venv/bin/activate
```
When the virtual environment is active, the command line should show: 
```
(.venv) user/path_to_file
```
Then, you will need to install the packages necessary to run the code within the virtual environment:
```
.venv/bin/python -m pip install numpy pandas requests biopython pymemesuite
```
Then, clone the repository into your using the following command in the terminal:
```
git clone https://github.com/annalundeen/CompBioProject
```
Then navigate into the project directory: 
```
cd CompBioProject/
```
Next, drag the newly created virtual environment into the project directory. Navigate to your home directory, locate the subdirectory ".venv", and drag it into "CompBioProject".

Once your virtual environment is set up, packages are installed, and you are in the correct directory, you can run the comprehensive python script: 
```
.venv/bin/python combined_motif_analysis.py
```
