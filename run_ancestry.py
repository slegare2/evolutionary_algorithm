#! /usr/bin/python3

import sys
sys.path.append('/home/slegare/IRIF/storydepth/evolutionary_algorithm')
from ancestry import RuleAncestry

past_dir = "past_generations"  # Directory where to read past models.
ance_dir = "ancestry"          # Directory where to write ancestry results

ancestry = RuleAncestry("best/best_100-1799.ka", "A binds B 2",
                         past_dir, ance_dir)

