DATABASE_NAME = 'multicom_3_sequence.db'
PDB_NAME = "multicom"

PDB_NAME_FILE = PDB_NAME + ".pdb"


import sqlite3
conn = sqlite3.connect(DATABASE_NAME)

c = conn.cursor()

c.execute("DROP TABLE IF EXISTS sequences")
c.execute('''CREATE TABLE sequences (seq text,
              acid_a text, acid_b text, acid_c text,
              acid_a_phi real, acid_a_psi real,
              acid_b_phi real, acid_b_psi real,
              acid_c_phi real, acid_c_psi real)''')

import Bio.PDB
for model in Bio.PDB.PDBParser().get_structure(PDB_NAME, PDB_NAME_FILE) :
    for chain in model :
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        
        for poly_index, poly in enumerate(polypeptides):

            sequence = poly.get_sequence()
            phi_psi_data = poly.get_phi_psi_list()
            
            number_acids = len(phi_psi_data)
            sql_input_data = []
            
            i = 2 #throw away the first and last amino acid
            while i<number_acids-2:
                
                seq = sequence[i-1], sequence[i], sequence[i+1]
                three_seq = [seq[0]+seq[1]+seq[2]]
                angles = phi_psi_data[i-1], phi_psi_data[i], phi_psi_data[i+1]                
                angle_data = [angles[0][0], angles[0][1], angles[1][0], angles[1][1], angles[2][0], angles[2][1]]
                
                sql_input_data += [[three_seq[0], seq[0], seq[1], seq[2], angle_data[0], angle_data[1], angle_data[2], angle_data[3], angle_data[4], angle_data[5]]]
                
                i+=1

            #print sql_input_data

            c.executemany('INSERT INTO sequences VALUES (?, ?,?,?, ?,?, ?,?, ?,?)', sql_input_data)

conn.commit()
conn.close()