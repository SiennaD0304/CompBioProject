## Code Overview
The entire analysis pipeline for this project is written in the combined_motif_analysis.py, which inputs a list of proteins and UniProt IDs to run BLAST, multiple sequence alignment, and MEME Suite FIMO to generate candidate motif matches. These candidate motifs are then used to generate PyMOL scripts to visualize each motif for putative interaction with the scorpion toxin LqhIII. 

## Project Setup
First, clone this repository into your environment using the following command in the terminal:
```
git clone https://github.com/annalundeen/CompBioProject
```
Then navigate into the project directory: 
```
cd CompBioProject/
```

## Virtual Environment and Running Code
This code requires a virtual environment so that packages can be installed in the project folder instead of globally. Documentation on setting up and managing environments in VS code can be found here: 
https://code.visualstudio.com/docs/python/environments.

First, make sure that you are in your project directory. Then, simply use the following command in the terminal, where your desired environment name can take the place of .venv: 
```
python3 -m venv .venv
```
Once the virtual environment is set up in your project folder, make sure that it is activated using the command:
```
source .venv/bin/activate
```
When the virtual environment is active, the terminal prompt should show: 
```
(.venv) user/path_to_project_directory
```
Then, you will need to install the packages necessary to run the code within the virtual environment:
```
.venv/bin/python -m pip install numpy pandas requests biopython pymemesuite
```
Once your virtual environment is set up, packages are installed, and you are in the correct directory, you can run the comprehensive python script: 
```
.venv/bin/python combined_motif_analysis.py
```

## Inputs and Expected Outputs
