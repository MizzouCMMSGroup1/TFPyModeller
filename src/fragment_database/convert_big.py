
import os
import subprocess

import input_files as input_data

for i in input_data.INPUT_PDB:
    #./rama input_pdb/1A1D_A.pdb -o output_rama/1A1D_A.txt
    file = i.split('.')
    file_name = file[0]
    os.system("./rama pdb_out/" + file_name + ".pdb -o output_big/" + file_name + ".txt")
    
    
    