INPUT_SEQUENCE = 'MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV'

DATABASE_NAME = 'multicom_3_sequence.db'

PDB_OUT = "PDB_OUT.txt"


import sqlite3
conn = sqlite3.connect(DATABASE_NAME)
conn.row_factory = sqlite3.Row

cursor = conn.cursor()

file = open(PDB_OUT, "w")

def pdb_append(acid, phi, psi):
    file.write("{0} {1:5} {2:5}\n".format(acid, phi, psi))

def append_sequence(a,b,c):
    sequence = a+b+c
    #print sequence
    
    cursor.execute("SELECT * FROM sequences WHERE seq = '%s'" % sequence)
    data = cursor.fetchone()
    
    if data == None:
        return
        
    pdb_append(a,data["acid_a_phi"],data["acid_a_psi"])
    pdb_append(b,data["acid_b_phi"],data["acid_b_psi"])
    pdb_append(c,data["acid_c_phi"],data["acid_c_psi"])
    


sequence_length = len(INPUT_SEQUENCE)

i = 0
while i<sequence_length-2:

    append_sequence(INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2])    
    i+=3
    
if (sequence_length%3)!=0:
    #need some code to deal with last one/two acids
    print sequence_length%3, "chars left"
    
    
# print len(PDB_OUT)
# print PDB_OUT

file.close()
