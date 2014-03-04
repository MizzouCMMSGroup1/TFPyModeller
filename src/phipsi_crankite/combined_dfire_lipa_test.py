INPUT_SEQUENCE = 'MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV'

DATABASE_NAME = 'multicom_3_sequence.db'

PDB_PREPEND_NAME = "simulation_"

NUMBER_SIMULATIONS = 2


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

def append_sequence(pdb_number, a,b,c):
    sequence = a+b+c
    #print sequence
    
    cursor.execute("SELECT * FROM sequences WHERE seq = '%s'" % sequence)
    data = cursor.fetchone()
    
    if data == None:
        return
        
    pdb_append(pdb_number,a,data["acid_a_phi"],data["acid_a_psi"])
    pdb_append(pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
    pdb_append(pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    

def build_model(pdb_number):

    sequence_length = len(INPUT_SEQUENCE)

    i = 0
    while i<sequence_length-2:

        append_sequence(pdb_number, INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2])    
        i+=3
    
    if (sequence_length%3)!=0:
        #need some code to deal with last one/two acids
        print sequence_length%3, "chars left"

def score_model(pdb_number):
    PDB_OUT = phipsi_file_name(pdb_number)
    PDB_OUT_LIPA = pdb_lipa_name(pdb_number)
    
    lipa_convert = os.system("./lipa " + PDB_OUT + " -o " + PDB_OUT_LIPA)
    #dfire_score = os.system("./dDFIRE " + PDB_OUT_LIPA)

    output = subprocess.check_output("./dDFIRE " + PDB_OUT_LIPA, shell=True)
    
    print "dfire_score: ", output



def main():

    for i in range(1,NUMBER_SIMULATIONS+1):
        build_model(i)
        
    for i in range(1,NUMBER_SIMULATIONS+1):
        score_model(i)


if __name__ == "__main__":
    main()