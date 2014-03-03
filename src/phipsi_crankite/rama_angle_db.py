RAMA_NAME = "multicom_rama.txt"
DATABASE_NAME = 'multicom_3_sequence.db'

import sqlite3
conn = sqlite3.connect(DATABASE_NAME)

c = conn.cursor()

c.execute("DROP TABLE IF EXISTS sequences")
c.execute('''CREATE TABLE sequences (seq text,
              acid_a text, acid_b text, acid_c text,
              acid_a_phi real, acid_a_psi real,
              acid_b_phi real, acid_b_psi real,
              acid_c_phi real, acid_c_psi real)''')

sequence_data = []
sql_input_data = []

with open(RAMA_NAME) as f:
	sequence_data = f.readlines()

sequence = sequence_data[1:] #first line is column data

number_acids = len(sequence_data)

i = 2
while i<number_acids-3:
    seq = sequence[i-1][0], sequence[i][0], sequence[i+1][0]
    three_seq = [seq[0]+seq[1]+seq[2]]
    angles_phi = float(sequence[i-1][2:8]), float(sequence[i][2:8]), float(sequence[i+1][2:8])
    angles_psi = float(sequence[i-1][9:15]), float(sequence[i][9:15]), float(sequence[i+1][9:15])
    angle_data = [angles_phi[0], angles_psi[0], angles_phi[1], angles_psi[1], angles_phi[2], angles_psi[2]]
    
    sql_input_data += [[three_seq[0], seq[0], seq[1], seq[2], angle_data[0], angle_data[1], angle_data[2], angle_data[3], angle_data[4], angle_data[5]]]
    
    #print three_seq, angles_phi, angles_psi
    i += 1

c.executemany('INSERT INTO sequences VALUES (?, ?,?,?, ?,?, ?,?, ?,?)', sql_input_data)

conn.commit()
conn.close()
