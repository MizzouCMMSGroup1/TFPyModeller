'''current test sequence'''
INPUT_SEQUENCE = 'MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV'

'''fragment database'''
DATABASE_NAME = 'fragment_database/300PDB_nine_sequences.db'
DATABASE_NAME = 'fragment_database/17K_PDB_nine_sequences.db'

'''store working files in temp_out'''
TMP_DIR = "temp_out/"

'''store working files in temp_out subdirectory/start of each file'''
PDB_PREPEND_NAME = TMP_DIR+"simulation_"

'''number of hill climbing simulations to run'''
NUMBER_SIMULATIONS = 1000

#some imports we need
import os
import subprocess
import random
import math
import copy

#blosum62
from blosum62 import Blosum62
blosum = Blosum62()
bmatrix = blosum.matrix

#keep track of our db
import sqlite3
conn = sqlite3.connect(DATABASE_NAME)
conn.row_factory = sqlite3.Row

cursor = conn.cursor()

#utility functions
def phipsi_file_name(phipsi_number):
    '''return phipsi/lipa file name for input number'''
    return PDB_PREPEND_NAME + str(phipsi_number) +".txt"
    
def pdb_lipa_name(pdb_number):
    '''return PDB file name for input number'''
    return PDB_PREPEND_NAME + "lipa_" + str(pdb_number) +".pdb"

def stats_save_name():
    return TMP_DIR+"stats" +".txt"

def stats_append(sequence_count,trial_number,score):
    '''add data to stats file'''
    
    STATS_OUT = stats_save_name()
    
    file = open(STATS_OUT, "a")
    file.write("{0}, {1}, {2}\n".format(sequence_count, trial_number, score))
    file.close()

#code for initialization of protein randomly from sequence database
def phipsi_append(phipsi_number, acid, phi, psi):
    '''add data to phipsi file'''
    
    PHIPSI_OUT = phipsi_file_name(phipsi_number)
    
    file = open(PHIPSI_OUT, "a")
    file.write("{0} {1:5} {2:5}\n".format(acid, phi, psi))
    file.close()

def append_sequence(pdb_number, a,b,c, d,e,f, g,h,i, skip_sequences=0):
    '''find sequence from database for each chunk of three residues, add to phipsi file'''
    
    sequence = a+b+c+d+e+f+g+h+i
    #print sequence
    
    cursor.execute("SELECT * FROM sequences WHERE seq = '%s' ORDER BY RANDOM() LIMIT 1;" % sequence)
    data = cursor.fetchone()
    
    if data == None:
        #need a better LIKE clause here (BLOSUM scores?)
        cursor.execute("SELECT * FROM sequences WHERE seq LIKE '%s_______%s' ORDER BY RANDOM() LIMIT 1;" % (a, i))
        data = cursor.fetchone()
        
        if data == None:
            "WARNING: no replacement found!"
            exit()
            return
        
        print "couldn't find sequence:", sequence, "using replacement:", data["seq"]
    
    #using modular arithmetic to only update the last few residues when needed
    
    if skip_sequences<1:
        phipsi_append(pdb_number,a,data["acid_a_phi"],data["acid_a_psi"])
    if skip_sequences<2:
        phipsi_append(pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
    if skip_sequences<3:
        phipsi_append(pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    if skip_sequences<4:
        phipsi_append(pdb_number,d,data["acid_d_phi"],data["acid_d_psi"])
    if skip_sequences<5:
        phipsi_append(pdb_number,e,data["acid_e_phi"],data["acid_e_psi"])
    if skip_sequences<6:
        phipsi_append(pdb_number,f,data["acid_f_phi"],data["acid_f_psi"])
    if skip_sequences<7:
        phipsi_append(pdb_number,g,data["acid_g_phi"],data["acid_g_psi"])
    if skip_sequences<8:
        phipsi_append(pdb_number,h,data["acid_h_phi"],data["acid_h_psi"])
    if skip_sequences<9:
        phipsi_append(pdb_number,i,data["acid_i_phi"],data["acid_i_psi"])

    

def build_model(pdb_number):
    '''build initial model for each chunk of three residues'''

    sequence_length = len(INPUT_SEQUENCE)

    i = 0
    while i<sequence_length-8:
        append_sequence(pdb_number, INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2],
                                    INPUT_SEQUENCE[i+3], INPUT_SEQUENCE[i+4], INPUT_SEQUENCE[i+5],
                                    INPUT_SEQUENCE[i+6], INPUT_SEQUENCE[i+7], INPUT_SEQUENCE[i+8] )    
        i+=9
    
    #return
    #some code to deal with last one/two acids
    if (sequence_length%9>0):
        append_sequence(pdb_number, INPUT_SEQUENCE[sequence_length-9], INPUT_SEQUENCE[sequence_length-8], INPUT_SEQUENCE[sequence_length-7],
                                    INPUT_SEQUENCE[sequence_length-6], INPUT_SEQUENCE[sequence_length-5], INPUT_SEQUENCE[sequence_length-4],
                                    INPUT_SEQUENCE[sequence_length-3], INPUT_SEQUENCE[sequence_length-2], INPUT_SEQUENCE[sequence_length-1],
                                    9-sequence_length%9)


#code for randomly replacing protein sub-sequences (from sequence database)
def phipsi_replace(old_model_number, new_model_number, a, a_phi, a_psi, b, b_phi, b_psi, c, c_phi, c_psi,
                                                       d, d_phi, d_psi, e, e_phi, e_psi, f, f_phi, f_psi,
                                                       g, g_phi, g_psi, h, h_phi, h_psi, i, i_phi, i_psi,
                                                       sequence_offset=-1):
    '''copy over old phipsi file to new one and add new sequence data'''
    if (sequence_offset<0):
        print "bad things happened"
        exit()
    else:
        print "seq offset:", sequence_offset, "old model:", old_model_number, "new model:", new_model_number

    PHIPSI_IN = phipsi_file_name(old_model_number)
    PHIPSI_OUT = phipsi_file_name(new_model_number)
    
    lines = [line.strip() for line in open(PHIPSI_IN)]
    
    a_out = "{0} {1:5} {2:5}\n".format(a, a_phi, a_psi)
    b_out = "{0} {1:5} {2:5}\n".format(b, b_phi, b_psi)
    c_out = "{0} {1:5} {2:5}\n".format(c, c_phi, c_psi)
    d_out = "{0} {1:5} {2:5}\n".format(d, d_phi, d_psi)
    e_out = "{0} {1:5} {2:5}\n".format(e, e_phi, e_psi)
    f_out = "{0} {1:5} {2:5}\n".format(f, f_phi, f_psi)
    g_out = "{0} {1:5} {2:5}\n".format(g, g_phi, g_psi)
    h_out = "{0} {1:5} {2:5}\n".format(h, h_phi, h_psi)
    i_out = "{0} {1:5} {2:5}\n".format(i, i_phi, i_psi)
    
    out_file = open(PHIPSI_OUT, "a")
    i = 0
    numlines = len(lines)
    
    while i < sequence_offset:
        out_file.write(lines[i]+"\n")
        i+=1
    
    out_file.write(a_out)
    out_file.write(b_out)
    out_file.write(c_out)
    out_file.write(d_out)
    out_file.write(e_out)
    out_file.write(f_out)
    out_file.write(g_out)
    out_file.write(h_out)
    out_file.write(i_out)
    i+=9
    
    while i < numlines:
        out_file.write(lines[i]+"\n")
        i+=1
    
    out_file.close()


def replace_sequence(old_model_number, new_model_number, a,b,c, d,e,f, g,h,i, sequence_offset=-1):
    '''find a new replacement sequence for an input'''
    if (sequence_offset<0):
        print "bad things happened"
        exit()
    
    #sequence = a+b+c+d+e+f+g+h+i
    #print sequence
    
    base_sequence = [a,b,c,d,e,f,g,h,i]
    sequence = copy.copy(base_sequence)
    
    #print sequence
    cursor.execute("SELECT * FROM sequences WHERE seq = '%s' ORDER BY RANDOM() LIMIT 1;" % ''.join(sequence))
    data = cursor.fetchone()
    
    j = 0
    k = 0
    
    while data == None:
        query = ''.join(sequence)
        cursor.execute("SELECT * FROM sequences WHERE seq = '%s' ORDER BY RANDOM() LIMIT 1;" % (query))
        data = cursor.fetchone()
        if data == None:
          sequence = copy.copy(base_sequence)
          sequence[j%9] = bmatrix[base_sequence[j%9]][0][0]
          if j%9 == 0 and j > 0:
            base_sequence[k%9] = bmatrix[base_sequence[k%9]][0][0]
            k = k+1
          print("no match found. using",''.join(sequence),"instead")
          j += 1
        if j > 17: # after three cycles we're too far away from the original sequence
            #need a better LIKE clause here (BLOSUM scores?)
            cursor.execute("SELECT * FROM sequences WHERE seq LIKE '%s_______%s' ORDER BY RANDOM() LIMIT 1;" % (a, i))
            data = cursor.fetchone()
    
            if data == None:
                "WARNING: no replacement found!"
                exit()
    
            print "couldn't find sequence:", sequence, "using replacement:", data["seq"]
    
    #print "skipping search"
    #return
    phipsi_replace(old_model_number, new_model_number,a,data["acid_a_phi"],data["acid_a_psi"],
                                                      b,data["acid_b_phi"],data["acid_b_psi"],
                                                      c,data["acid_c_phi"],data["acid_c_psi"],
                                                      d,data["acid_d_phi"],data["acid_d_psi"],
                                                      e,data["acid_e_phi"],data["acid_e_psi"],
                                                      f,data["acid_f_phi"],data["acid_f_psi"],
                                                      g,data["acid_g_phi"],data["acid_g_psi"],
                                                      h,data["acid_h_phi"],data["acid_h_psi"],
                                                      i,data["acid_i_phi"],data["acid_i_psi"],
                                                       sequence_offset)


def randomize_model(old_model_number, new_model_number):
    '''pick three residues at random and replace them'''

    sequence_length = len(INPUT_SEQUENCE)
    
    random_offset = random.randint(0,sequence_length-9) #can't select beyond last residue

    replace_sequence(old_model_number, new_model_number, INPUT_SEQUENCE[random_offset], INPUT_SEQUENCE[random_offset+1], INPUT_SEQUENCE[random_offset+2],
                                                         INPUT_SEQUENCE[random_offset+3], INPUT_SEQUENCE[random_offset+4], INPUT_SEQUENCE[random_offset+5],
                                                         INPUT_SEQUENCE[random_offset+6], INPUT_SEQUENCE[random_offset+7], INPUT_SEQUENCE[random_offset+8],
                                                         random_offset)

#global utility function to build/score pdb file
def build_score_pdb_model(pdb_number):
    '''send phipsi to lipa to generated pdb, score with dfire, return score (most negative == best)'''
    PDB_OUT = phipsi_file_name(pdb_number)
    PDB_OUT_LIPA = pdb_lipa_name(pdb_number)
    
    lipa_convert = os.system("./lipa " + PDB_OUT + " -o " + PDB_OUT_LIPA)

    dfire_output = subprocess.check_output("./dDFIRE " + PDB_OUT_LIPA, shell=True)

    dfire_output_split = dfire_output.split(":")
    dfire_output_split_nums = dfire_output_split[1].split(" ")

    print "pdb:", PDB_OUT_LIPA, "dfire_score:", float(dfire_output_split_nums[1])

    return float(dfire_output_split_nums[1])


#temperature calculator. non-linear decrease
def sigmoid_temperature(k):
  return -5000/(1 + math.exp(-k/200)) + 5000

def linear_temperature(k):
  return (-2500/1000)*k + 2500

#init and run hill climbing, print out results
def main():

    #cleanup
    if not os.path.exists(TMP_DIR):
        print("Making temp directory")
        os.mkdir(TMP_DIR)
    else:
        print("Clearing temp directory")
        for root,dirs,files in os.walk(TMP_DIR,topdown=False):
            for name in files:
                os.remove(os.path.join(root,name))
    '''
    #cleanup
    for i in range(1,NUMBER_SIMULATIONS+1):
        if os.path.exists(phipsi_file_name(i)): os.rm()
        os.system("rm " + phipsi_file_name(i))
        os.system("rm " + pdb_lipa_name(i))
    '''

    #build initial random model
    build_model(0)
    
    best_dfire_score = 1e6 #best is lowest
    best_dfire_number = -1

    old_model_number = 0

    #test models
    #model 0 is the initial sequence
    for i in range(1,NUMBER_SIMULATIONS+1,1):
        
        #T = sigmoid_temperature(i)
        T = linear_temperature(i)
        
        randomize_model(old_model_number, i)
        
        i_score = build_score_pdb_model(i)
        
        stats_append(9,i,i_score)
        
        score_diff = i_score - best_dfire_score
        
        if score_diff > 0:
            #print("score diff", score_diff)
            score_diff = score_diff if score_diff < 5*T else 5*T
            prob_to_accept = math.exp(-100*score_diff/T)
            print("probability to accept:", prob_to_accept)
            if prob_to_accept < random.random():
                #print("probability not enough")
                continue
            print("accepting some randomness")

        #print("scoring")
        #if i_score < best_dfire_score:
        best_dfire_score = i_score
        best_dfire_number = i
        old_model_number = i
    
    
    #save best
    stats_append(9000,best_dfire_number,best_dfire_score)
    
    print "best model found:", best_dfire_number, "score:", best_dfire_score


if __name__ == "__main__":
    main()

