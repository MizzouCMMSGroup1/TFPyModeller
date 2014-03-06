DATABASE_NAME = '17K_PDB_nine_sequences.db'

import input_files as input_data
#test
#RAMA_FILES = ['1A1L_A.txt']

import sqlite3
conn = sqlite3.connect(DATABASE_NAME)

c = conn.cursor()

c.execute("DROP TABLE IF EXISTS sequences")
c.execute('''CREATE TABLE sequences (seq text,
              acid_a_phi real, acid_a_psi real,
              acid_b_phi real, acid_b_psi real,
              acid_c_phi real, acid_c_psi real,
              acid_d_phi real, acid_d_psi real,
              acid_e_phi real, acid_e_psi real,
              acid_f_phi real, acid_f_psi real,
              acid_g_phi real, acid_g_psi real,
              acid_h_phi real, acid_h_psi real,
              acid_i_phi real, acid_i_psi real  )''')

def parse_rama_file(RAMA_NAME):
    
    sequence_data = []
    sql_input_data = []

    with open('output_big/' + RAMA_NAME) as f:
    	sequence_data = f.readlines()

    sequence = sequence_data[1:] #first line is column data

    number_acids = len(sequence)

    nine_seq = []

    i = 1 #throw away the first and last amino acid
    while i < number_acids-9:
        nine_seq = []
        angles_phi = []
        angles_psi = []
        for s in range(0,9):
            nine_seq += sequence[i+s][0]
            angles_phi += [float(sequence[i+s][2:8])]
            angles_psi += [float(sequence[i+s][9:15])]
            
        #print nine_seq
        #print angles_phi, angles_psi
        
        sql_input_data += [[''.join(nine_seq),
                angles_phi[0], angles_psi[0], 
                angles_phi[1], angles_psi[1], 
                angles_phi[2], angles_psi[2], 
                angles_phi[3], angles_psi[3], 
                angles_phi[4], angles_psi[4], 
                angles_phi[5], angles_psi[5], 
                angles_phi[6], angles_psi[6], 
                angles_phi[7], angles_psi[7], 
                angles_phi[8], angles_psi[8]  ]]
        
        i += 1

    c.executemany('INSERT INTO sequences VALUES (?, ?,?, ?,?, ?,?, ?,?, ?,?, ?,?,  ?,?, ?,?, ?,?)', sql_input_data)


def main():

    for f in input_data.INPUT_PDB:
        file = f.split('.')
        file_name = file[0]
        parse_rama_file(file_name+".txt")
        
    conn.commit()
    conn.close()


if __name__ == "__main__":
    main()
