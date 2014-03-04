INPUT_SEQUENCE = 'MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV'

DATABASE_NAME = 'fragment_database/300PDB_three_sequences.db'

PDB_PREPEND_NAME = "temp_out/simulation_"

NUMBER_SIMULATIONS = 5


import os
import subprocess

import sqlite3
conn = sqlite3.connect(DATABASE_NAME)
conn.row_factory = sqlite3.Row

cursor = conn.cursor()

def phipsi_file_name(pdb_number):
    return PDB_PREPEND_NAME + str(pdb_number) +".txt"
    
def pdb_lipa_name(pdb_number):
    return PDB_PREPEND_NAME + "lipa_" + str(pdb_number) +".pdb"

def pdb_append(pdb_number, acid, phi, psi):
    PDB_OUT = phipsi_file_name(pdb_number)
    
    file = open(PDB_OUT, "a")
    file.write("{0} {1:5} {2:5}\n".format(acid, phi, psi))
    file.close()

def append_sequence(pdb_number, a,b,c, skip_sequences=0):
    sequence = a+b+c
    #print sequence
    
    cursor.execute("SELECT * FROM sequences WHERE seq = '%s' ORDER BY RANDOM() LIMIT 1;" % sequence)
    data = cursor.fetchone()
    
    if data == None:
        #need a better LIKE clause here (BLOSUM scores?)
        cursor.execute("SELECT * FROM sequences WHERE seq LIKE '%s_%s' ORDER BY RANDOM() LIMIT 1;" % (a, c))
        data = cursor.fetchone()
        
        if data == None:
            "WARNING: no replacement found!"
            exit()
            return
        
        print "skipping sequence:", sequence, "replacement:", data["seq"]
    
    #using modular arithmetic to only update the last 1 or 2 residues when needed
    if skip_sequences==0:
        pdb_append(pdb_number,a,data["acid_a_phi"],data["acid_a_psi"])
        pdb_append(pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
        pdb_append(pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    
    elif skip_sequences==2:
        pdb_append(pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
        pdb_append(pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    
    elif skip_sequences==1:
        pdb_append(pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    

def build_model(pdb_number):

    sequence_length = len(INPUT_SEQUENCE)

    i = 0
    while i<sequence_length-2:
        append_sequence(pdb_number, INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2])    
        i+=3
    
    #some code to deal with last one/two acids
    if (sequence_length%3>0):
        append_sequence(pdb_number, INPUT_SEQUENCE[sequence_length-3], INPUT_SEQUENCE[sequence_length-2], INPUT_SEQUENCE[sequence_length-1], sequence_length%3)
    

def score_model(pdb_number):
    PDB_OUT = phipsi_file_name(pdb_number)
    PDB_OUT_LIPA = pdb_lipa_name(pdb_number)
    
    lipa_convert = os.system("./lipa " + PDB_OUT + " -o " + PDB_OUT_LIPA)

    dfire_output = subprocess.check_output("./dDFIRE " + PDB_OUT_LIPA, shell=True)

    dfire_output_split = dfire_output.split(":")
    dfire_output_split_nums = dfire_output_split[1].split(" ")

    print "pdb:", PDB_OUT_LIPA, "dfire_score:", float(dfire_output_split_nums[1])

    return float(dfire_output_split_nums[1])

def main():

    #cleanup
    for i in range(1,NUMBER_SIMULATIONS+1):
        os.system("rm " + phipsi_file_name(i))
        os.system("rm " + pdb_lipa_name(i))

    #build models
    for i in range(1,NUMBER_SIMULATIONS+1):
        build_model(i)
    
    best_dfire_score = 1e6 #best is lowest
    best_dfire_number = -1

    #test models
    for i in range(1,NUMBER_SIMULATIONS+1):
        i_score = score_model(i)
        if i_score < best_dfire_score:
            best_dfire_score = i_score
            best_dfire_number = i
            
    print "best model found:", best_dfire_number, "score:", best_dfire_score


if __name__ == "__main__":
    main()

