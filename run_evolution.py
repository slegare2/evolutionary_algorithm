#! /usr/bin/python3

import sys
sys.path.append('/home/slegare/IRIF/storydepth/evolutionary_algorithm')
from evolutionary import Evolutionary
from mutaterule import MutateRule, MutateRate


# ==== Set Evolutionary Algorithm Parameters. ====

# Simulation length params.
num_gens   = 100   # Number of generations to simulate.
sim_time   = 10    # Simulation time of each simulation.
out_period = 0.01  # Time period between writings in the cvs output file.
                   # Should be smaller or equal to time_ave_frac/sim_time such
                   # that enough points are used to compute average amounts.

# Fitness and selection params.
selection_type = "fit"  # Can be "fit" (for fitness proportional) or "rank".
elite_frac     = 0.1    # Fraction of the population considered as elite.
penalty_weight = 0.2    # Weight of the fitness penalty for number of rules.
time_ave_frac  = 0.1    # Final portion of time series used to compute
                        # observable averages.

# Mutation params.
mutation_prob = 1     # Probability that a child is mutated. Problem right now if not 1: clones erase their father and population size decreases.
refine_prob   = 0.1   # Probability that a mutation is a refinement.
                      # The mutation is otherwise a rate change.
max_bonds     = 3     # Maximum number of bonds allowed from primary binding
                      # partner in refinements

# Available rates params.
binary_rates = [0, 0.1]
#binary_rates = [0, 0.00001, 0.001, 0.1]
unary_rates = [0, 1000]

# Input directories.
start_dir = "starting_pop"     # Directory containinf the initial models.
pop_dir = "current_generation" # Directory containing the models to run.

# Output directories.
next_dir = "next_generation"   # Directory where to write the next generation.
past_dir = "past_generations"  # Directory where to copy past models.
out_dir  = "time_series"       # Directory where to write output time series.
fit_dir  = "fitness"           # Directory where to write computed fitness.
best_dir = "best"              # Directory where to keep the best individual
                               # of each generation
story_dir = "stories"          # Directory where to write stories.

# Interaction matrix file.
matrix_file = "../inter_matrix_full6.txt"

# ================================================

# ---------- Run Evolutionary Algorithm ----------

evolution = Evolutionary(num_gens, sim_time, out_period,
                         selection_type, elite_frac,
                         penalty_weight, time_ave_frac,
                         mutation_prob, refine_prob, max_bonds,
                         binary_rates, unary_rates,
                         start_dir, pop_dir, next_dir, past_dir, out_dir,
                         fit_dir, best_dir, story_dir, matrix_file)
evolution.evolve()
#evolution.compute_stories([0, 20, 100])

# ------------------------------------------------
