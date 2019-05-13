#! /usr/bin/python3

""" 
Creates an initial population of different individuals with mutated rules.
"""

import os
import sys
import shutil
import random
import filecmp
sys.path.append('/home/slegare/IRIF/storydepth/evolutionary_algorithm')
from mutaterule import MutateRule, MutateRate


# User parameters.
max_bonds = 3   # Maximum number of edges from the binding partner.
num_indiv = 20  # Number of individuals in the initial population.
binary_rates = [0, 0.1]
unary_rates = [0, 1000]
refine_prob = 0.2

# Interaction matrix file and Kappa model.
matrix_file = "inter_matrix_full6.txt"
kappa_file = "abc-0.ka"

# Some manipulations on the Kappa file.
dash = kappa_file.index("-")
prefix = kappa_file[:dash+1]
shutil.copyfile(kappa_file, "init_pop/{}".format(kappa_file))

for model_num in range(1, num_indiv):
    is_different = False
    found_identical = False
    while is_different == False:
        refine_val = random.uniform(0, 1)
        if refine_val <= refine_prob: # Refine a Rule.
            mutated_model_tmp = MutateRule(kappa_file, matrix_file, max_bonds)
            tmp_out = "init-tmp.ka"
            output_content = open(tmp_out, "w")
            output_content.write(mutated_model_tmp.new_file)
            output_content.close()
            mut_lines = mutated_model_tmp.mutated_lines
            mutated_model = MutateRate(tmp_out, binary_rates, unary_rates,
                                       select_lines=mut_lines,
                                       unchanged_father=True)
        else: # Change a rate.
            mutated_model = MutateRate(kappa_file, binary_rates, unary_rates)
        # Write new model.
        output_path = "init_pop/{}{}.ka".format(prefix, model_num)
        output_content = open(output_path, "w")
        output_content.write(mutated_model.new_file)
        output_content.close()
        if refine_val <= refine_prob:
            os.remove(tmp_out)
        # Check if the new model is not identical to one that already exists.
        for previous_num in range(0, model_num):
            is_different = True
            previous_path = "init_pop/{}{}.ka".format(prefix, previous_num)
            previous_content = open(previous_path, "r").readlines()
            output_content = open(output_path, "r").readlines()
            if len(previous_content) == len(output_content):
                is_different = False
                for i in range(len(output_content)):
                    previous_line = previous_content[i]
                    output_line = output_content[i]
                    if previous_line != output_line:
                        is_different = True
                        break
            if is_different == False:
                print("File {}{}.ka is identical to {}{}.ka. Making a new model."
                      .format(prefix, model_num, prefix, previous_num))
                found_identical = True
                break
        if is_different == True:
            if found_identical == True:
                now = " now"
            else:
                now = ""
            print("File {}{}.ka is{} different from all previous models."
                  .format(prefix, model_num, now))
 
