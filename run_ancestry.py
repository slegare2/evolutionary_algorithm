#! /usr/bin/python3

import sys
sys.path.append('/home/slegare/IRIF/storydepth/evolutionary_algorithm')
from ancestry import Ancestry

fit_dir  = "fitness"           # Directory where to read models and fitness.
past_dir = "past_generations"  # Directory where to read past models.
ance_dir = "ancestry"          # Directory where to write ancestry results.

ancestry = Ancestry(fit_dir, past_dir, ance_dir)

