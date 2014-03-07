'''current test sequence'''
INPUT_SEQUENCE = 'MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPSQAIQDALKKLYPNATAIKWEQKGVYYVADCQADGREKEVWFDANANWLMTETELNSINNLPPAVLTAFMESSYNNWVVDDVVILEYPNEPSTEFVVTVEQGKKVDLYFSEGGGLLHEKDVTNGDDTHWPRV'

INPUT_SEQUENCE = 'MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPP'

'''fragment database'''
THREE_SEQUENCE_DATABASE_NAME = 'fragment_database/17K_PDB_three_sequences.db'
NINE_SEQUENCE_DATABASE_NAME = 'fragment_database/17K_PDB_nine_sequences.db'

'''directory to store all target files'''
TARGET_DIR = "targets/"

'''name for files'''
TARGET_NAME = "T9999"

'''number of hill climbing simulations to run'''
NUMBER_SIMULATIONS = 1000

#some imports we need
import os
import subprocess
import random
import math
import copy
import shutil

#blosum62
from blosum62 import Blosum62
blosum = Blosum62()
bmatrix = blosum.matrix

#keep track of our db
import sqlite3
conn3 = sqlite3.connect(THREE_SEQUENCE_DATABASE_NAME)
conn3.row_factory = sqlite3.Row
conn9 = sqlite3.connect(NINE_SEQUENCE_DATABASE_NAME)
conn9.row_factory = sqlite3.Row

cursor3 = conn3.cursor()
cursor9 = conn9.cursor()

#utility functions
def phipsi_file_name(seq_length,phipsi_number):
    '''return phipsi/lipa file name for input number'''
    return TARGET_DIR+TARGET_NAME+"/phipsi/"+ "sequence%d-%d" % (seq_length,phipsi_number) +".txt"
    
def pdb_lipa_name(seq_length,pdb_number):
    '''return PDB file name for input number'''
    return TARGET_DIR+TARGET_NAME+"/lipa/"+ "sequence%d-%d" % (seq_length,pdb_number) +".pdb"

def pdb_save_name(pdb_number):
    return TARGET_DIR+TARGET_NAME+"/pdbs/"+ "sequence%d" % (pdb_number) +".pdb"

#code for initialization of protein randomly from sequence database
def phipsi_append(seq_length,phipsi_number, acid, phi, psi):
    '''add data to phipsi file'''
    
    PHIPSI_OUT = phipsi_file_name(seq_length,phipsi_number)
    
    file = open(PHIPSI_OUT, "a")
    file.write("{0} {1:>6} {2:>6}\n".format(acid, phi, psi))
    file.close()

'''
Replace9 Functions
'''

def append_sequence9(pdb_number, a,b,c, d,e,f, g,h,i, skip_sequences=0, head=False, tail=False):
    '''find sequence from database for each chunk of three residues, add to phipsi file'''
    
    base_sequence = [a,b,c,d,e,f,g,h,i]
    sequence = copy.copy(base_sequence)
    #print sequence
    print("finding sequence matching:", ''.join(sequence))

    data = None
    j = 0
    k = 0

    while data == None:
        query = ''.join(sequence)
        cursor9.execute("SELECT * FROM sequences WHERE seq = '%s' ORDER BY RANDOM() LIMIT 1;" % (query))
        data = cursor9.fetchone()
        if data == None:
          sequence = copy.copy(base_sequence)
          sequence[j%9] = bmatrix[base_sequence[j%9]][0][0]
          if j%9 == 0 and j > 0:
            base_sequence[k%9] = bmatrix[base_sequence[k%9]][0][0]
            k = k+1
          print("no match found. using",''.join(sequence),"instead")
          j += 1
        if j > 17: # after three cycles we're too far away from the original sequence
          if skip_sequences > 0:
                if skip_sequences < 3:
                    append_sequence3(pdb_number, a,b,c, skip_sequences=skip_sequences)
                    append_sequence3(pdb_number, d,e,f)
                    append_sequence3(pdb_number, g,h,i, tail=tail)
                elif skip_sequences < 6:
                    append_sequence3(pdb_number, d,e,f, skip_sequences=skip_sequences-3)
                    append_sequence3(pdb_number, g,h,i, tail=tail)
                elif skip_sequences < 9:
                    append_sequence3(pdb_number, g,h,i, skip_sequences=skip_sequences-6, tail=tail)
          else:
              append_sequence3(pdb_number, a,b,c, head=head)
              append_sequence3(pdb_number, d,e,f)
              append_sequence3(pdb_number, g,h,i)
          return

    #using modular arithmetic to only update the last few residues when needed

    if skip_sequences<1:
        if head:
            phipsi_append(9,pdb_number,'nan',data["acid_a_psi"])
        else:
            phipsi_append(9,pdb_number,a,data["acid_a_phi"],data["acid_a_psi"])
    if skip_sequences<2:
        phipsi_append(9,pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
    if skip_sequences<3:
        phipsi_append(9,pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    if skip_sequences<4:
        phipsi_append(9,pdb_number,d,data["acid_d_phi"],data["acid_d_psi"])
    if skip_sequences<5:
        phipsi_append(9,pdb_number,e,data["acid_e_phi"],data["acid_e_psi"])
    if skip_sequences<6:
        phipsi_append(9,pdb_number,f,data["acid_f_phi"],data["acid_f_psi"])
    if skip_sequences<7:
        phipsi_append(9,pdb_number,g,data["acid_g_phi"],data["acid_g_psi"])
    if skip_sequences<8:
        phipsi_append(9,pdb_number,h,data["acid_h_phi"],data["acid_h_psi"])
    if skip_sequences<9:
        if tail:
            phipsi_append(9,pdb_number,i,data["acid_i_phi"],'nan')
        else:
            phipsi_append(9,pdb_number,i,data["acid_i_phi"],data["acid_i_psi"])

    

def build_model9(pdb_number):
    '''build initial model for each chunk of three residues'''

    sequence_length = len(INPUT_SEQUENCE)

    i = 0
    append_sequence9(pdb_number, INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2],
                                 INPUT_SEQUENCE[i+3], INPUT_SEQUENCE[i+4], INPUT_SEQUENCE[i+5],
                                 INPUT_SEQUENCE[i+6], INPUT_SEQUENCE[i+7], INPUT_SEQUENCE[i+8],
                                 head=True)

    i = 9
    while i<sequence_length-9:
        append_sequence9(pdb_number, INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2],
                                     INPUT_SEQUENCE[i+3], INPUT_SEQUENCE[i+4], INPUT_SEQUENCE[i+5],
                                     INPUT_SEQUENCE[i+6], INPUT_SEQUENCE[i+7], INPUT_SEQUENCE[i+8] )    
        i+=9
    
    #return
    #some code to deal with last one/two acids
    if (sequence_length%9>0):
        append_sequence9(pdb_number, INPUT_SEQUENCE[sequence_length-9], INPUT_SEQUENCE[sequence_length-8], INPUT_SEQUENCE[sequence_length-7],
                                     INPUT_SEQUENCE[sequence_length-6], INPUT_SEQUENCE[sequence_length-5], INPUT_SEQUENCE[sequence_length-4],
                                     INPUT_SEQUENCE[sequence_length-3], INPUT_SEQUENCE[sequence_length-2], INPUT_SEQUENCE[sequence_length-1],
                                     skip_sequences=sequence_length%9,tail=True)

    src_file = phipsi_file_name(3,pdb_number)
    dest_file = phipsi_file_name(9,pdb_number)
    shutil.copy(src_file,dest_file)
    os.remove(src_file)

#code for randomly replacing protein sub-sequences (from sequence database)
def phipsi_replace9(old_model_number, new_model_number, a, a_phi, a_psi, b, b_phi, b_psi, c, c_phi, c_psi,
                                                       d, d_phi, d_psi, e, e_phi, e_psi, f, f_phi, f_psi,
                                                       g, g_phi, g_psi, h, h_phi, h_psi, i, i_phi, i_psi,
                                                       sequence_offset=-1):
    '''copy over old phipsi file to new one and add new sequence data'''
    if (sequence_offset<0):
        raise BaseException("phipsi_replace3 requires a sequence offset >= 0.", sequence_offset, "was given")
        exit()
    else:
        print("seq offset:", sequence_offset, "old model:", old_model_number, "new model:", new_model_number)

    PHIPSI_IN = phipsi_file_name(9,old_model_number)
    # clone off the file in case of updating itself
    shutil.copy(PHIPSI,'temp.txt')
    PHIPSI_IN = 'temp.txt'
    PHIPSI_OUT = phipsi_file_name(9,new_model_number)
    
    lines = [line.strip() for line in open(PHIPSI_IN)]
    numlines = len(lines)
    
    if seq_offset == 0:
        a_out = "{0} {1:>6} {2:>6}\n".format(a, 'nan', a_psi)
    else:
        a_out = "{0} {1:>6} {2:>6}\n".format(a, a_phi, a_psi)

    b_out = "{0} {1:>6} {2:>6}\n".format(b, b_phi, b_psi)
    c_out = "{0} {1:>6} {2:>6}\n".format(c, c_phi, c_psi)
    d_out = "{0} {1:>6} {2:>6}\n".format(d, d_phi, d_psi)
    e_out = "{0} {1:>6} {2:>6}\n".format(e, e_phi, e_psi)
    f_out = "{0} {1:>6} {2:>6}\n".format(f, f_phi, f_psi)
    g_out = "{0} {1:>6} {2:>6}\n".format(g, g_phi, g_psi)
    h_out = "{0} {1:>6} {2:>6}\n".format(h, h_phi, h_psi)

    if seq_offset+8 == numlines-1:
        i_out = "{0} {1:>6} {2:>6}\n".format(i, i_phi, 'nan')
    else:
        i_out = "{0} {1:>6} {2:>6}\n".format(i, i_phi, i_psi)
    
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
    
    os.remove(PHIPSI_IN)


def replace_sequence9(old_model_number, new_model_number, a,b,c, d,e,f, g,h,i, sequence_offset=-1):
    '''find a new replacement sequence for an input'''
    if (sequence_offset<0):
        raise BaseException("replace_sequence9 requires a sequence offset >= 0.", sequence_offset, "was given")
        exit()
    
    base_sequence = [a,b,c,d,e,f,g,h,i]
    sequence = copy.copy(base_sequence)
    #print sequence
    print("finding sequence matching:", ''.join(sequence))

    data = None
    j = 0
    k = 0

    while data == None:
        query = ''.join(sequence)
        cursor9.execute("SELECT * FROM sequences WHERE seq = '%s' ORDER BY RANDOM() LIMIT 1;" % (query))
        data = cursor9.fetchone()
        if data == None:
          sequence = copy.copy(base_sequence)
          sequence[j%9] = bmatrix[base_sequence[j%9]][0][0]
          if j%9 == 0 and j > 0:
            base_sequence[k%9] = bmatrix[base_sequence[k%9]][0][0]
            k = k+1
          print("no match found. using",''.join(sequence),"instead")
          j += 1
        if j > 17: # after three cycles we're too far away from the original sequence
          print("reverting to sequence length 3")
          # Create a length3 copy of the file in order
          # to use the replace_sequence3 method, which
          # only takes a length3 phipsi file
          src_file = phipsi_file_name(9,old_model_number)
          dest_file = phipsi_file_name(3,old_model_number)
          shutil.copy(src_file,dest_file)
          replace_sequence3(old_model_number, new_model_number, a,b,c, sequence_offset)
          replace_sequence3(new_model_number, new_model_number, d,e,f, sequence_offset+3)
          replace_sequence3(new_model_number, new_model_number, g,h,i, sequence_offset+6)
          # Swap the file back and cleanup your mess
          os.remove(dest_file)
          src_file = phipsi_file_name(3,new_model_number)
          dest_file = phipsi_file_name(9,new_model_number)
          shutil.copy(src_file,dest_file)
          os.remove(src_file)
          return

    phipsi_replace9(old_model_number, new_model_number,a,data["acid_a_phi"],data["acid_a_psi"],
                                                      b,data["acid_b_phi"],data["acid_b_psi"],
                                                      c,data["acid_c_phi"],data["acid_c_psi"],
                                                      d,data["acid_d_phi"],data["acid_d_psi"],
                                                      e,data["acid_e_phi"],data["acid_e_psi"],
                                                      f,data["acid_f_phi"],data["acid_f_psi"],
                                                      g,data["acid_g_phi"],data["acid_g_psi"],
                                                      h,data["acid_h_phi"],data["acid_h_psi"],
                                                      i,data["acid_i_phi"],data["acid_i_psi"],
                                                       sequence_offset)


def randomize_model9(old_model_number, new_model_number):
    '''pick three residues at random and replace them'''

    sequence_length = len(INPUT_SEQUENCE)
    
    random_offset = random.randint(0,sequence_length-9) #can't select beyond last residue

    replace_sequence9(old_model_number, new_model_number, INPUT_SEQUENCE[random_offset], INPUT_SEQUENCE[random_offset+1], INPUT_SEQUENCE[random_offset+2],
                                                         INPUT_SEQUENCE[random_offset+3], INPUT_SEQUENCE[random_offset+4], INPUT_SEQUENCE[random_offset+5],
                                                         INPUT_SEQUENCE[random_offset+6], INPUT_SEQUENCE[random_offset+7], INPUT_SEQUENCE[random_offset+8],
                                                         random_offset)

'''
Replace3 Functions
'''

def append_sequence3(pdb_number, a,b,c, skip_sequences=0,head=False,tail=False):
    '''find sequence from database for each chunk of three residues, add to phipsi file'''
    
    sequence = [a,b,c]
    #print sequence
    print("finding sequence matching:", ''.join(sequence))
    
    data = None
    i = 0
    
    while data == None:
      query = ''.join(sequence)
      cursor3.execute("SELECT * FROM sequences WHERE seq= '%s' ORDER BY RANDOM() LIMIT 1;" % query)
      data = cursor3.fetchone()
      sequence[i%3] = bmatrix[sequence[i%3]][0][0]
      i+=1
    
    #using modular arithmetic to only update the last 1 or 2 residues when needed
    if skip_sequences==0:
        if head:
          phipsi_append(3,pdb_number,a,'nan',data["acid_a_psi"])
        else:
          phipsi_append(3,pdb_number,a,data["acid_a_phi"],data["acid_a_psi"])
          
        phipsi_append(3,pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
        
        if tail:
          phipsi_append(3,pdb_number,c,data["acid_c_phi"],'nan')
        else:
          phipsi_append(3,pdb_number,c,data["acid_c_phi"],data["acid_c_psi"])
    
    elif skip_sequences==2:
        phipsi_append(3,pdb_number,b,data["acid_b_phi"],data["acid_b_psi"])
        phipsi_append(pdb_number,c,data["acid_c_phi"],'nan')
    
    elif skip_sequences==1:
        phipsi_append(3,pdb_number,c,data["acid_c_phi"],'nan')
    

def build_model3(pdb_number):
    '''build initial model for each chunk of three residues'''

    sequence_length = len(INPUT_SEQUENCE)

    append_sequence3(pdb_number, INPUT_SEQUENCE[0], INPUT_SEQUENCE[1], INPUT_SEQUENCE[2],head=True)
    i=3
    while i<sequence_length-3:
        append_sequence3(pdb_number, INPUT_SEQUENCE[i], INPUT_SEQUENCE[i+1], INPUT_SEQUENCE[i+2])    
        i+=3
    
    #some code to deal with last one/two acids
    append_sequence3(pdb_number, INPUT_SEQUENCE[sequence_length-3], INPUT_SEQUENCE[sequence_length-2], INPUT_SEQUENCE[sequence_length-1], skip_sequences=sequence_length%3,tail=True)


#code for randomly replacing protein sub-sequences (from sequence database)
def phipsi_replace3(old_model_number, new_model_number, a, a_phi, a_psi, b, b_phi, b_psi, c, c_phi, c_psi, sequence_offset=-1):
    '''copy over old phipsi file to new one and add new sequence data'''
    if (sequence_offset<0):
        raise BaseException("phipsi_replace3 requires a sequence offset >= 0.", sequence_offset, "was given")
        exit()
    else:
        print("seq offset:", sequence_offset, "old model:", old_model_number, "new model:", new_model_number)

    PHIPSI_IN = phipsi_file_name(3,old_model_number)
    # clone off the file in case of updating itself
    shutil.copy(PHIPSI_IN,'temp.txt')
    PHIPSI_IN = 'temp.txt'
    PHIPSI_OUT = phipsi_file_name(3,new_model_number)
    
    lines = [line.strip() for line in open(PHIPSI_IN)]
    numlines = len(lines)
    
    if sequence_offset == 0:
      a_out =  "{0} {1:>6} {2:>6}\n".format(a, 'nan', a_psi)
    else:
      a_out = "{0} {1:>6} {2:>6}\n".format(a, a_phi, a_psi)
    
    b_out = "{0} {1:>6} {2:>6}\n".format(b, b_phi, b_psi)
    
    if sequence_offset+2 == numlines-1:
      c_out = "{0} {1:>6} {2:>6}\n".format(c, c_phi, 'nan')
    else:
      c_out = "{0} {1:>6} {2:>6}\n".format(c, c_phi, c_psi)
    
    out_file = open(PHIPSI_OUT, "a")
    i = 0
    
    while i < sequence_offset:
        out_file.write(lines[i]+"\n")
        i+=1
    
    out_file.write(a_out)
    out_file.write(b_out)
    out_file.write(c_out)
    i+=3
    
    while i < numlines:
        out_file.write(lines[i]+"\n")
        i+=1
    
    out_file.close()
    
    os.remove(PHIPSI_IN)


def replace_sequence3(old_model_number, new_model_number, a,b,c, sequence_offset=-1):
    '''find a new replacement sequence for an input'''
    if (sequence_offset<0):
        raise BaseException("replace_sequence3 requires a sequence offset >= 0.", sequence_offset, "was given")
        exit()
    
    sequence = [a,b,c]
    #print sequence
    print("finding sequence matching:", ''.join(sequence))
    
    data = None
    i = 0
    
    while data == None:
      query = ''.join(sequence)
      cursor3.execute("SELECT * FROM sequences WHERE seq= '%s' ORDER BY RANDOM() LIMIT 1;" % query)
      data = cursor3.fetchone()
      sequence[i%3] = bmatrix[sequence[i%3]][0][0]
      i+=1
    
    phipsi_replace3(old_model_number, new_model_number,a,data["acid_a_phi"],data["acid_a_psi"],b,data["acid_b_phi"],data["acid_b_psi"],c,data["acid_c_phi"],data["acid_c_psi"], sequence_offset)


def randomize_model3(old_model_number, new_model_number):
    '''pick three residues at random and replace them'''

    sequence_length = len(INPUT_SEQUENCE)
    
    random_offset = random.randint(0,sequence_length-3) #can't select beyond last residue

    replace_sequence3(old_model_number, new_model_number, INPUT_SEQUENCE[random_offset], INPUT_SEQUENCE[random_offset+1], INPUT_SEQUENCE[random_offset+2], random_offset)

#global utility function to build/score pdb file
def build_score_pdb_model(seq_length,pdb_number):
    '''send phipsi to lipa to generated pdb, score with dfire, return score (most negative == best)'''
    PDB_OUT = phipsi_file_name(seq_length,pdb_number)
    PDB_OUT_LIPA = pdb_lipa_name(seq_length,pdb_number)
    
    lipa_convert = os.system("./lipa " + PDB_OUT + " -o " + PDB_OUT_LIPA)
    
    score = -1
    
    while score < 0:
      dfire_output = subprocess.check_output("./dDFIRE " + PDB_OUT_LIPA, shell=True)
  
      dfire_output_split = dfire_output.split(":")
      dfire_output_split_nums = dfire_output_split[1].split(" ")
  
      score = float(dfire_output_split_nums[1])

    #print("pdb:", PDB_OUT_LIPA, "dfire_score:", float(dfire_output_split_nums[1])
    print("pdb:", pdb_number, "score:", float(score))

    return score

#temperature calculator. non-linear decrease
def sigmoid_temperature(k):
  return -5000/(1 + math.exp(-k/200)) + 5000

def linear_temperature(k):
  return (-2500/1000)*k + 2500
  
#init and run hill climbing, print out results
def main(temperature=sigmoid_temperature):

    #cleanup
    if not os.path.exists(TARGET_DIR+TARGET_NAME):
        print("Making temp directory")
        os.mkdir(TARGET_DIR+TARGET_NAME)
    else:
        print("Clearing temp directory")
        for root,dirs,files in os.walk(TARGET_DIR+TARGET_NAME,topdown=False):
            for name in dirs:
                shutil.rmtree(os.path.join(root,name))

    os.mkdir(TARGET_DIR+TARGET_NAME+"/phipsi")
    os.mkdir(TARGET_DIR+TARGET_NAME+"/lipa")
    os.mkdir(TARGET_DIR+TARGET_NAME+"/pdbs")

    #build initial random model
    build_model9(0)
    
    best_dfire_score = 1e6 #best is lowest
    best_dfire_number = -1

    old_model_number = 0
    num_pdbs = 1

    best_model = old_model_number
    best_score = best_dfire_score
    best_dfire = best_dfire_number

    for k in range(1,NUMBER_SIMULATIONS+1):
        # Grab temperature
        T = temperature(k)
        # Create neighbor
        randomize_model9(best_model, k)
        # Calculate score. Lower is better
        neighbor_score = build_score_pdb_model(9,k)
        # Calculate score difference
        score_diff = neighbor_score - best_score
        # If score_diff is above 0, chance it
        if score_diff > 0:
            # Grab the probability of accepting a bad neighbor
            prob_to_accept = math.exp(-100*score_diff/T)
            print("probability to accept:", prob_to_accept)
            # Reject if roll is above probability
            if prob_to_accept < random.random():
                continue
            print("accepted neighbor anyway")
        best_model = k
        best_dfire = k
        best_score = neighbor_score

        src_name = pdb_lipa_name(9,best_model)
        dest_name = pdb_save_name(num_pdbs)
        shutil.copy(src_name,dest_name)
        num_pdbs += 1

    for k in range(1,NUMBER_SIMULATIONS+1):
        # Grab temperature
        T = temperature(k)
        # Create neighbor
        randomize_model3(best_model, k)
        # Calculate score. Lower is better
        neighbor_score = build_score_pdb_model(3,k)
        # Calculate score difference
        score_diff = neighbor_score - best_score
        # If score_diff is above 0, chance it
        if score_diff > 0:
            # Grab the probability of accepting a bad neighbor
            prob_to_accept = math.exp(-100*score_diff/T)
            print("probability to accept:", prob_to_accept)
            # Reject if roll is above probability
            if prob_to_accept < random.random():
                continue
            print("accepted neighbor anyway")
        best_model = k
        best_dfire = k
        best_score = neighbor_score

        src_name = pdb_lipa_name(3,best_model)
        dest_name = pdb_save_name(num_pdbs)
        shutil.copy(src_name,dest_name)
        num_pdbs += 1

    print("best model found:", best_model, "score:", best_score)

def run(target,sequence=INPUT_SEQUENCE,simulations=1000,temperature='sigmoid'):
  global TARGET_NAME
  global INPUT_SEQUENCE
  global NUMBER_SIMULATIONS

  print("Runnig Simulated Annleaing for target:",target)
  if target:
    TARGET_NAME = target

  if sequence:
    INPUT_SEQUENCE = sequence

  if simulations:
    NUMBER_SIMULATION = simulations

  if temperature in ['sigmoid','linear']:
    if temperature == 'sigmoide':
      main(sigmoid_temperature)
    if temperature == 'linear':
      main(linear_temperature)
  else:
    main(sigmoid_temperature)

def runargs(args):
    target = TARGET_NAME
    sequence = INPUT_SEQUENCE
    simulations = NUMBER_SIMULATIONS
    temperature = 'sigmoid'
    if args[0]:
        target = args[0]
    if args[1]:
        sequence = args[1]
    if args[2]:
        simulations = args[2]
    if args[3]:
        temperature = args[3]
    run(target,sequence,simulations,temperature)
    
if __name__ == "__main__":
    main()

