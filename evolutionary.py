#! /usr/bin/python3

""" Evolutionary algorithm. """

import os
import time
import shutil
import random
import string
import math
import json
import kappy
import graphviz
from mutaterule import MutateRule, MutateRate


class Evolutionary:
    """ Evolutionary algorithm. """

    def __init__(self, num_gens, sim_time, replicates, out_period,
                 selection_type, elite_frac,
                 penalty_weight, time_ave_frac, score_operation,
                 mutation_prob, refine_prob, max_bonds,
                 binary_rates, unary_rates,
                 start_dir, pop_dir, next_dir, past_dir, out_dir,
                 fit_dir, best_dir, story_dir, din_dir, matrix_file):
        """ Initialize Evolutionary class. """

        self.num_gens = num_gens
        self.sim_time = sim_time
        self.replicates = replicates
        self.out_period = float(out_period)
        self.selection_type = selection_type
        self.elite_frac = elite_frac
        self.penalty_weight = penalty_weight
        self.time_ave_frac = time_ave_frac
        self.score_operation = score_operation
        self.mutation_prob = mutation_prob
        self.refine_prob = refine_prob
        self.max_bonds = max_bonds
        self.binary_rates = binary_rates
        self.unary_rates = unary_rates
        self.start_dir = start_dir
        self.pop_dir = pop_dir
        self.next_dir = next_dir
        self.past_dir = past_dir
        self.out_dir = out_dir
        self.fit_dir = fit_dir
        self.best_dir = best_dir
        self.story_dir = story_dir
        self.din_dir = din_dir
        self.matrix_file = matrix_file
        self.letters = list(string.ascii_lowercase)
        # Some checks.
        if self.selection_type != "fit" and self.selection_type != "rank":
            raise Exception("Parameter selection_type should be either "
                            "'fit' or 'rank'.")


    def evolve(self):
        """
        Loop over the evolutionary algorithm to
        produce num_gens generations.
        """

        self.track_generation()
        while self.generation < self.num_gens:
            self.generation += 1
            self.model_list = self.get_files(self.pop_dir)
            self.check_survivors()
            self.simulation()
            self.compute_fitness()
            self.selection()
            self.mutation()
            self.keep_best()
            self.replace_population()
            print("Generation {} complete.".format(self.generation))
            self.update_gen_num()
        self.summarize()


    def track_generation(self):
        """ Keep track of the last generation produced. """

        try:
            gen_file = open("generation_num.txt", "r").readlines()
            self.generation = int(gen_file[0])
        except:
            self.generation = 0


    def update_gen_num(self):
        """ Write the latest generation number to file generation_num.txt. """

        gen_file = open("generation_num.txt", "w")
        gen_file.write("{}".format(self.generation))
        gen_file.close()


    def get_files(self, directory):
        """
        Get all files from a given directory and sort them.
        Also cut the extention of the file.
        """

        tmp_list = os.listdir(directory)
        file_dicts = []
        for f in tmp_list:
            dash = f.rfind("-")
            dot = f.rfind(".")
            name = f[:dot]
            extension = f[dot+1:]
            if extension == "dat" and self.replicates > 1:
                number = int(f[dash+1:dot-1])
            else:
                number = int(f[dash+1:dot])
            file_dicts.append({"file": name, "num": number})
        sorted_dicts = sorted(file_dicts, key=lambda x: x["num"])
        file_list = []
        for d in sorted_dicts:
            file_list.append(d["file"])
        return file_list


    def check_survivors(self):
        """
        Get the list of models that survived unmutated (elites) from
        the previous generation.
        """

        self.survivors = []
        self.survivors_fit = []
        self.surv_num_rules = {}
        self.surv_time_averages = {}
        if self.generation > 1:
            file_path = "{}/fitness_gen-{}.txt".format(self.fit_dir,
                                                       self.generation-1)
            fitness_file = open(file_path).readlines()
            for i in range(2, len(fitness_file)):
                tokens = fitness_file[i].split()
                if tokens[0] in self.model_list:
                    surv_model = tokens[0]
                    surv_score = float(tokens[-3])
                    surv_pen = float(tokens[-2])
                    surv_fit = float(tokens[-1])
                    self.survivors.append(surv_model)
                    self.survivors_fit.append([surv_model, surv_score,
                                               surv_pen, surv_fit])
                    self.surv_num_rules[surv_model] = int(tokens[1])
                    obs_list = []
                    for i in range(2, len(tokens)-3):
                        obs_list.append(float(tokens[i]))
                    self.surv_time_averages[surv_model] = obs_list
        self.new_model_list = []
        for model in self.model_list:
            if model not in self.survivors:
                self.new_model_list.append(model)

    # ============== Simulation Section ===================

    def simulation(self):
        """ Run a KaSim simulation on each model found in pop_dir.  """

        self.clear_previous()
        for surv in self.survivors:
            input_path = "{}/{}.ka".format(self.pop_dir, surv)
            print("Skipping {} which was already simulated."
                  .format(input_path))
        for model in self.new_model_list:
            input_path = "{}/{}.ka".format(self.pop_dir, model)
            kappa_file = open(input_path, "r").read()
            for sim_num in range(self.replicates):
                if self.replicates > 1:
                    letter_id = self.letters[sim_num]
                    replicate = "{}{}".format(model, letter_id)
                    print("Running simulation {} on {}."
                          .format(sim_num+1, input_path))
                else:
                    replicate = model
                    print("Running simulation on {}."
                          .format(input_path))
                client = kappy.KappaStd()
                client.add_model_string(kappa_file)
                client.project_parse()
                sim_params = kappy.SimulationParameter(
                    pause_condition="[T] > {}".format(self.sim_time),
                    plot_period=self.out_period)
                client.simulation_start(sim_params)
                while client.get_is_sim_running(): time.sleep(0.1)
                results = client.simulation_plot()
                # Write simulation results to output file.
                self.write_output(replicate, results)
                client.shutdown()


    def clear_previous(self):
        """
        Clear time series from the simulations of the previous
        generation.
        """

        output_list = self.get_files(self.out_dir)
        for output in output_list:
            file_path = "{}/{}.dat".format(self.out_dir, output)
            #print("removing file {}".format(file_path))
            os.remove(file_path)


    def write_output(self, model, results):
        """ Write time series to output directory. """

        output_path = "{}/{}.dat".format(self.out_dir, model)
        output_file = open(output_path, "w")
        # Write header.
        header = results["legend"]
        for title in header:
            output_file.write("{}   ".format(title))
        output_file.write("\n")
        # Write time series.
        time_series = results["series"]
        n = len(time_series)
        for i in range(n-1, -1, -1):
            entry = time_series[i]
            for value in entry:
                output_file.write("{:10.4f}   ".format(value))
            output_file.write("\n")

    # ============= Simulation Section End ================

    # ----------- Fitness Calculation Section -------------

    def compute_fitness(self):
        """ Compute fitness of all simulated model. """

        self.get_obs_weights()
        self.get_num_rules()
        self.compute_time_averages()
        self.fitness = self.survivors_fit.copy()
        for model in self.new_model_list:
            # Compute additive score.
            add_score = 0
            for i in range(self.num_obs):
                w = self.weights[i]
                n = self.time_averages[model][i]
                add_score += w*n
            # Compute multiplicative score.
            mult_score = 1
            for i in range(self.num_obs):
                w = self.weights[i]
                n = self.time_averages[model][i]
                mult_score *= w*n
            if self.score_operation == "add":
                score = add_score
            elif self.score_operation == "mult":
                score = mult_score
            elif self.score_operation == "addmult":
                score = add_score + mult_score
            else:
                raise Exception('score_operation should be either "add", '
                                '"mult" or "addmult"')
            # Compute penalty.
            if self.penalty_weight != 0.0:
                p = self.penalty_weight
                n = self.num_rules[model] - self.min_num_rules
                penalty = p*n
            else:
                penalty = 0
            # Compute fitness
            fitness = score - penalty
            data = [model, score, penalty, fitness]
            self.fitness.append(data)
        self.rank_models()
        self.write_fitness()


    def get_obs_weights(self):
        """
        Get observables and their weights from first model.
        Weights should be the same for all models.
        """

        input_path = "{}/{}.ka".format(self.pop_dir, self.model_list[0])
        kappa_file = open(input_path, "r").readlines()
        self.observables = []
        self.weights = []
        for line in kappa_file:
            if line[:5] == "%obs:":
                begin_name = line.index("'")
                end_name = line.rfind("'")
                name = line[begin_name+1:end_name]
                tokens = line.split()
                weight = float(tokens[-1])
                self.observables.append(name)
                self.weights.append(weight)
        self.num_obs = len(self.observables)


    def get_num_rules(self):
        """ Get the number of rules for all models. """

        self.num_rules = self.surv_num_rules.copy()
        for model in self.new_model_list:
            input_path = "{}/{}.ka".format(self.pop_dir, model)
            kappa_file = open(input_path, "r").readlines()
            n_rules = 0
            counting = False
            for line in kappa_file:
                if "// Unary binding rules."in line:
                    break
                if counting == True:
                    n_rules += 1
                if "// Binary binding rules." in line:
                    counting = True
            self.num_rules[model] = n_rules
        # Get the number of rules of the minimal model.
        if self.penalty_weight != 0.0:
            starting_files = self.get_files(self.start_dir)
            for f in starting_files:
                dash = f.rfind("-")
                model_num = f[dash+1:]
                if model_num == "0":
                    self.model_0_file = f
                    break
            input_path = "{}/{}.ka".format(self.start_dir, self.model_0_file)
            kappa_file = open(input_path, "r").readlines()
            n_rules = 0
            counting = False
            for line in kappa_file:
                if "// Unary binding rules."in line:
                    break
                if counting == True:
                    n_rules += 1
                if "// Binary binding rules." in line:
                    counting = True
            self.min_num_rules = n_rules
            


    def compute_time_averages(self):
        """
        Compute the time average of each observable
        for each simulation.
        """
        
        self.time_averages = self.surv_time_averages.copy()
        for model in self.new_model_list:
            replicate_averages = []
            # Compute observable averages of each simulation.
            for sim_num in range(self.replicates):
                if self.replicates > 1:
                    letter_id = self.letters[sim_num]
                    replicate = "{}{}".format(model, letter_id)
                else:
                    replicate = model
                input_path = "{}/{}.dat".format(self.out_dir, replicate)
                time_series = open(input_path, "r").readlines()
                n_plot = len(time_series)
                n_entries = n_plot - 1
                n_read = int(n_entries * self.time_ave_frac)
                read_start = n_plot - n_read
                obs_averages = []
                for i in range(1, self.num_obs+1):
                    summ = 0.0
                    for j in range(read_start, n_plot):
                        tokens = time_series[j].split()
                        summ += float(tokens[i])
                    if n_read != 0:
                        ave = summ / n_read
                    else:
                        ave = 0
                    obs_averages.append(ave)
                replicate_averages.append(obs_averages)
            # Compute the global averages.
            glob_averages = []
            for i in range(self.num_obs):
                summ = 0.0
                for replica in replicate_averages:
                    summ += replica[i]
                ave = summ / self.replicates
                glob_averages.append(ave)
            self.time_averages[model] = glob_averages


    def rank_models(self):
        """ Rank each simulated model based on fitness. """

        self.ranked_models = sorted(self.fitness, key=lambda x: x[-1],
                                    reverse=True)
        # Shuffle models with same fitness.
        group_list = [ [self.ranked_models[0]] ]
        group_index = 0
        for i in range(1, len(self.ranked_models)):
            prev_fit = self.ranked_models[i-1][-1]
            curr_fit = self.ranked_models[i][-1]
            if curr_fit == prev_fit:
                group_list[group_index].append(self.ranked_models[i])
            else:
                group_list.append([self.ranked_models[i]])
                group_index += 1
        self.shuffled_models = []
        for group in group_list:
            random.shuffle(group)
            for model in group:
                self.shuffled_models.append(model)


    def write_fitness(self):
        """ Write fitness file for given generation. """

        output_path = "{}/fitness_gen-{}.txt".format(self.fit_dir,
                                                     self.generation)
        fitness_file = open(output_path, "w")
        # Write header.
        fitness_file.write("{:25} ".format("Model"))
        fitness_file.write("Rules   ")
        for obs in self.observables:
            fitness_file.write("{:15}".format(obs))
        fitness_file.write("  Score       Penalty     Fitness")
        fitness_file.write("\n\n")
        # Write data.
        for data in self.ranked_models:
            model = data[0]
            n_rules = self.num_rules[model]
            fitness_file.write("{:25} {:3}".format(model, n_rules))
            for obs in self.time_averages[model]:
                fitness_file.write("     {:7}   ".format(obs))
            fitness_file.write("  ")
            for entry in data[1:]:
                fitness_file.write("{:12.3f}".format(entry))
            fitness_file.write("\n")

    # --------- Fitness Calculation Section End -----------

    # ~~~~~~~~~~~~~ Selection Section ~~~~~~~~~~~~~~~~~~~~~

    def selection(self):
        """ Select individuals to be the parents of the next generation. """

        self.pop_size = len(self.model_list)
        self.n_elite = int(self.pop_size*self.elite_frac)
        self.n_fill = self.pop_size - self.n_elite
        self.elitism() # Does not do anything if elite_frac = 0
        if self.selection_type == "fit":
            self.fitness_proportional()
        elif self.selection_type == "rank":
            self.rank_selection()


    def elitism(self):
        """
        Pick the top elite_frac portion of the current population as elites.
        Elite individuals are directly copied into the next generation
        without mutation, for the moment.
        A future version will allow mutations but replace back children with
        parents if children have a lower fitness than their parents.
        """

        for i in range(self.n_elite):
            model = self.shuffled_models[i][0]
            input_path = "{}/{}.ka".format(self.pop_dir, model)
            output_path = "{}/{}.ka".format(self.next_dir, model)
            shutil.copyfile(input_path, output_path)

    
    def fitness_proportional(self):
        """
        Fitness proportional (roulette wheel) selection of individuals.
        Since negative fitness is possible, fitness is scaled by substracting
        the lowest fitness of the population. The lowest fitness individual
        hence has a zero probability of being selected, which is fine.
        """

        ordered_models = []
        ordered_fitness = []
        for entry in self.ranked_models:
            ordered_models.append(entry[0])
            ordered_fitness.append(entry[-1])
        min_fitness = min(ordered_fitness)
        scaled_fitness = []
        for fit in ordered_fitness:
            scaled_value = fit - min_fitness
            scaled_fitness.append(scaled_value)
        self.selected_models = random.choices(population=ordered_models,
                                              weights=scaled_fitness,
                                              k=self.n_fill)


    def rank_selection(self):
        """ Rank based selection of individuals. """
        
        ordered_models = []
        ranks = []
        rank = self.pop_size
        for entry in self.ranked_models:
            ordered_models.append(entry[0])
            ranks.append(rank)
            rank += -1
        self.selected_models = random.choices(population=ordered_models,
                                              weights=ranks,
                                              k=self.n_fill)

    # ~~~~~~~~~~~~~ Selection Section End ~~~~~~~~~~~~~~~~~ 

    # ............... Mutation Section ....................

    def mutation(self):
        """
        Mutate the selected individuals with probabilities given by
        mutation_prob and refine_prob.
        """

        self.track_model_num()
        for model in self.selected_models:
            input_path = "{}/{}.ka".format(self.pop_dir, model)
            mutate_val = random.uniform(0, 1)
            if mutate_val <= self.mutation_prob: # If true, mutate child.
                self.model_num += 1
                is_different = False
                found_identical = False
                while is_different == False:
                    refine_val = random.uniform(0, 1)
                    dash = model.rfind("-")
                    prefix = model[:dash+1]
                    if refine_val <= self.refine_prob: # Refine a Rule.
                        mutated_model_tmp = MutateRule(input_path,
                                                       self.matrix_file,
                                                       self.max_bonds)
                        tmp_out = "{}/{}tmp.ka".format(self.next_dir, prefix,
                                                       self.model_num)
                        output_content = open(tmp_out, "w")
                        output_content.write(mutated_model_tmp.new_file)
                        output_content.close()
                        mut_lines = mutated_model_tmp.mutated_lines
                        mutated_model = MutateRate(tmp_out, self.binary_rates,
                                                   self.unary_rates,
                                                   select_lines=mut_lines,
                                                   unchanged_father=True)
                    else: # Change a rate.
                        mutated_model = MutateRate(input_path, self.binary_rates,
                                                   self.unary_rates)
                    # Write mutated model to new file.
                    output_path = "{}/{}{}.ka".format(self.next_dir, prefix,
                                                      self.model_num)
                    output_content = open(output_path, "w")
                    output_content.write(mutated_model.new_file)
                    output_content.close()
                    # Clean temporary file.
                    if refine_val <= self.refine_prob:
                        os.remove(tmp_out)
                    # Check if the new model is not identical to one that
                    # already exists in next_dir.
                    created_files = self.get_files(self.next_dir)
                    created_files.remove("{}{}".format(prefix, self.model_num))
                    for created_file in created_files:
                        is_different = True
                        previous_path = "{}/{}.ka".format(self.next_dir, created_file)
                        previous_content = open(previous_path, "r").readlines()
                        output_content = open(output_path, "r").readlines()
                        if len(previous_content) == len(output_content):
                            is_different = False
                            # I start the following loop at 1 to skip the
                            # first line which contains the father.
                            for i in range(1, len(output_content)):
                                previous_line = previous_content[i]
                                output_line = output_content[i]
                                if previous_line != output_line:
                                    is_different = True
                                    break
                        if is_different == False:
                            #print("File {}{}.ka is identical to {}.ka. Making a new model."
                            #      .format( prefix, self.model_num, created_file))
                            found_identical = True
                            break
                    #if is_different == True:
                    #    if found_identical == True:
                    #        now = " now"
                    #    else:
                    #        now = ""
                    #    print("File {}{}.ka is{} different from all previous models."
                    #          .format(prefix, self.model_num, now))
            else: # If no mutation, just directly copy into next generation.
                output_path = "{}/{}.ka".format(self.next_dir, model)
                shutil.copyfile(input_path, output_path)
        self.update_model_num()


    def track_model_num(self):
        """ Keep track of the highest model number used. """

        if self.generation == 1:
            self.model_num = len(self.model_list)-1
        else:
            track_file = open("model_num.txt", "r").readlines()
            self.model_num = int(track_file[0])


    def update_model_num(self):
        """ Write the highest model number used to file model_num.txt. """

        track_file = open("model_num.txt", "w")
        track_file.write("{}".format(self.model_num))
        track_file.close()

    # .............. Mutation Section End .................

    def keep_best(self):
        """
        Keep the kappa file and output of the best individual of each
        generation.
        """

        best_model = self.ranked_models[0][0]
        dash = best_model.rfind("-")
        model_number = best_model[dash+1:]
        input_path = "{}/{}.ka".format(self.pop_dir, best_model)
        output_path = "{}/best_{}-{}.ka".format(self.best_dir,
                                                self.generation,
                                                model_number)
        shutil.copyfile(input_path, output_path)
        #input_data_path = "{}/{}.dat".format(self.out_dir, best_model)
        #output_data_path = "{}/best_ori_{}-{}.dat".format(self.best_dir,
        #                                              self.generation,
        #                                              model_number)
        #shutil.copyfile(input_data_path, output_data_path)


    def replace_population(self):
        """ Replace current generation with next generation. """

        # Move previous models to past_dir.
        for model in self.model_list:
            from_path = "{}/{}.ka".format(self.pop_dir, model)
            dash = model.rfind("-")
            num = model[dash+1:]
            to_path = "{}/gen_{}-{}.ka".format(self.past_dir,
                                               self.generation, num)
            #print("moving file {} to {}".format(from_path, to_path))
            shutil.move(from_path, to_path)
        # Copy next generation models into current generation.
        self.next_list = self.get_files(self.next_dir)
        for model in self.next_list:
            from_path = "{}/{}.ka".format(self.next_dir, model)
            to_path = "{}/{}.ka".format(self.pop_dir, model)
            #print("copying file {} to {}".format(from_path, to_path))
            shutil.copyfile(from_path, to_path)
        # Empty next generation directory.
        for model in self.next_list:
            file_path = "{}/{}.ka".format(self.next_dir, model)
            #print("removing file {}".format(file_path))
            os.remove(file_path)


    def summarize(self):
        """ Create a file to show the evolution of fitness over generations. """

        self.fit_list = self.get_files(self.fit_dir)
        self.generation_averages = []
        for fit in self.fit_list:
            self.compute_gen_ave(fit)
        output_file = open("fitness_summary.dat", "w")
        output_file.write("Gen      Fave        Fmax        Fmin"
                           "        Rave        Rmax        Rmin\n")
        for i in range(len(self.generation_averages)):
            output_file.write("{:3d} ".format(i+1))
            row = self.generation_averages[i]
            for val in row:
                output_file.write(" {:11.6f}".format(val))
            output_file.write("\n")


    def compute_gen_ave(self, fit_file):
        """ Compute the average values from one fitness file. """

        input_path = "{}/{}.txt".format(self.fit_dir, fit_file)
        fitness_file = open(input_path, "r").readlines()
        n_col = len(fitness_file[0].split())
        n_row = len(fitness_file)
        gen_ave = []
        # Fitness.
        val_list = []
        summ = 0
        for i in range(2, n_row):
            row = fitness_file[i]
            tokens = row.split()
            val = float(tokens[-1])
            summ += val
            val_list.append(val)
        gen_ave.append(summ / (n_row - 1))
        gen_ave.append(max(val_list))
        gen_ave.append(min(val_list))
        # Number of rules.
        val_list = []
        summ = 0
        for i in range(2, n_row):
            row = fitness_file[i]
            tokens = row.split()
            val = float(tokens[1])
            summ += val
            val_list.append(val)
        gen_ave.append(summ / (n_row - 1))
        gen_ave.append(max(val_list))
        gen_ave.append(min(val_list))
        self.generation_averages.append(gen_ave)

    # +++++++++++++ Story Section ++++++++++++++

    def compute_stories(self, selected_generations):
        """
        Compute the stories of each observable for the best model of selected
        generation. Generation 0 means the initial unrefined model.
        """

        self.selected_gens = selected_generations
        self.model_list = self.get_files(self.pop_dir)
        self.get_sel_models()
        self.create_directories()
        self.run_traces()
        self.run_kastor()
        self.get_obs_weights()
        self.sort_stories()
        self.draw_stories()


    def get_sel_models(self):
        """ Find the best model of the selected generations. """

        best_files = self.get_files(self.best_dir)
        best_models = []
        for f in best_files:
            if "ori" not in f:
                best_models.append(f)
        self.selected_models = []
        for model in best_models:
            underscore = model.index("_")
            dash = model.index("-")
            gen_num = int(model[underscore+1:dash])
            if gen_num in self.selected_gens:
                if model not in self.selected_models:
                    self.selected_models.append(model)
        if 0 in self.selected_gens:
            starting_files = self.get_files(self.start_dir)
            for f in starting_files:
                dash = f.rfind("-")
                model_num = f[dash+1:]
                if model_num == "0":
                    self.model_0_file = f
                    break
            

    def create_directories(self):
        """
        Create a directory for each observable and a subdirectory for each
        generation whose story is computed, if they do not already exist.
        """

        # Model 0 directory.        
        if 0 in self.selected_gens:
            model_0_path = "{}/{}".format(self.story_dir, self.model_0_file)
            if not os.path.exists(model_0_path):
                os.mkdir(model_0_path)
            input_path = "{}/{}.ka".format(self.start_dir, self.model_0_file)
            output_path = "{}/{}.ka".format(model_0_path, self.model_0_file)
            if not os.path.exists(output_path):
                shutil.copyfile(input_path, output_path)
        # Selected model directories.
        for model in self.selected_models:
            model_path = "{}/{}".format(self.story_dir, model)
            if not os.path.exists(model_path):
                os.mkdir(model_path)
            input_path = "{}/{}.ka".format(self.best_dir, model)
            output_path = "{}/{}.ka".format(model_path, model)
            if not os.path.exists(output_path):
                shutil.copyfile(input_path, output_path)
        # Add model_0 after this step.
        self.selected_models.append(self.model_0_file)


    def run_traces(self):
        """ Resimulate the best models but storing the trace. """

        for model in self.selected_models:
            input_path = "{}/{}/{}.ka".format(self.story_dir, model, model)
            output_path = "{}/{}/{}.dat".format(self.story_dir, model, model)
            trace_path = "{}/{}/{}.json".format(self.story_dir, model, model)
            command_line = ("KaSim -mode batch --no-log --no-log -u t -p {} "
                            "-l {} -i {} -o {} -trace {}"
                            .format(self.out_period, self.sim_time, input_path,
                                    output_path, trace_path))
            os.system(command_line)


    def run_kastor(self):
        """ Get stories using KaStor. """

        for model in self.selected_models:
            input_path = "{}/{}/{}.json".format(self.story_dir, model, model)
            output_path = "{}/{}/story".format(self.story_dir, model)
            command_line = ("KaStor --weak {} -o {}"
                            .format(input_path, output_path))
            os.system(command_line)

    
    def determine_observables(self, model):
        """
        Determine which observable is the event of interest in each story.
        """

        input_path = "{}/{}".format(self.story_dir, model)
        file_list = os.listdir(input_path)
        dot_list = []
        for f in file_list:
            if f[-4:] == ".dot":
                dot_list.append(f)
        self.obs_map = {}
        for i in range(len(dot_list)):
            file_name = "{}/{}/{}".format(self.story_dir, model, dot_list[i])
            underscore = file_name.rfind("_")
            story_id = file_name[underscore+1:-4]
            file_content = open(file_name, "r").readlines()
            last_label = 0
            for j in range(len(file_content)):
                line = file_content[j]
                if "label=" in line:
                    last_label = j
            obs_line = file_content[last_label]
            open_quote = obs_line.index('"')
            close_quote = obs_line.rfind('"')
            obs = obs_line[open_quote+1:close_quote]
            #if obs not in self.observables:
            #    raise Exception("Event of interest does not correspond to "
            #                    "a defined observable.")
            if obs not in self.obs_map.keys():
                self.obs_map[obs] = [story_id]
            else:
                self.obs_map[obs].append(story_id)


    def sort_stories(self):
        """ Count how many times each story occured. """

        self.model_dict = {}
        for model in self.selected_models:
            self.determine_observables(model)
            summary_name = ("{}/{}/story_Weakly_Summary.dat"
                           .format(self.story_dir, model))
            outfile_name = ("{}/{}/story_sort_weak.txt"
                           .format(self.story_dir, model))
            summary_file = open(summary_name, "r").readlines()
            outfile = open(outfile_name, "w")
            story_count = {}
            depths = {}
            for i in range(1, len(summary_file)):
                tokens = summary_file[i].split()
                story_id = tokens[0]
                if story_id not in story_count.keys():
                    story_count[story_id] = 1
                else:
                    story_count[story_id] = story_count[story_id] + 1
                depth = tokens[3]
                if story_id not in depths.keys():
                    depths[story_id] = depth
            split_obs = {}
            obs_sum = {}
            for obs in self.obs_map.keys():
                story_ids = self.obs_map[obs]
                numbers_list = []
                n = 0
                for story_id in story_ids:
                    numbers_list.append([story_id, story_count[story_id]])
                    n += int(story_count[story_id])
                sorted_numbers = sorted(numbers_list, key=lambda x: x[1], reverse= True)
                split_obs[obs] = sorted_numbers
                obs_sum[obs] = n
            outfile.write('Rank  id     N    frac   depth\n\n')
            for obs in split_obs.keys():
                outfile.write("{}\n\n".format(obs))
                story_list = split_obs[obs]
                rank = 0
                n = obs_sum[obs]
                for story in story_list:
                    story_id = story[0]
                    num = story[1]
                    frac = float(num)/n
                    rank += 1
                    depth = depths[story_id]
                    outfile.write("{:3}  {:3}  {:4}    {:5.3f}   {:3}\n"
                                  .format(rank, int(story_id), num, frac,
                                          depth))
                outfile.write("\n\n")
            self.model_dict[model] = split_obs


    def draw_stories(self):
        """
        Draw every story as a png.
        Also rename them in a meaningful way.
        """

        for model in self.selected_models:
            dir_path = "{}/{}".format(self.story_dir, model)
            split_obs = self.model_dict[model]
            for obs in split_obs.keys():
                story_list = split_obs[obs]
                rank = 0
                for story in story_list:
                    story_id = story[0]
                    num = story[1]
                    rank += 1
                    rank_str = self.add_zeros(rank, len(story_list))
                    input_path = ("{}/story_Weakly_{}.dot"
                                  .format(dir_path, story_id))
                    output_path = ("{}/{}_w_r{}_n{}_{}.png"
                                   .format(dir_path, obs, rank_str,
                                           num, story_id))
                    command_line = ("/usr/bin/dot -Tpng {} > {}"
                                    .format(input_path, output_path))
                    os.system(command_line)


    def add_zeros(self, rank, max_rank):
        """
        Add zeros in front of integers for them to appear ordered by 'ls'.
        """
        max_len = len(str(max_rank))
        rank_len = len(str(rank))
        to_add = max_len - rank_len
        r_str = ""
        for i in range(to_add):
            r_str += "0"
        r_str = r_str + str(rank)
        return r_str

    # +++++++++++ Story Section End ++++++++++++

    # """"""""""""" DIN section """"""""""""""""

    def compute_din(self, flux_time_frac):
        """
        Compute and draw flux maps of the best model of each generation.
        Unbinding and unimolecular rules are not considered at the moment.
        """

        self.flux_time_frac = flux_time_frac
        self.get_best_models()
        self.create_flux_maps()
        self.draw_flux_maps()


    def get_best_models(self):
        """ Get the best models to simulate. """

        best_files = self.get_files(self.best_dir)
        self.best_dict = {}
        for f in best_files:
            dash = f.rfind("-")
            underscore = f.rfind("_")
            model_num = f[dash+1:]
            gen_num = f[underscore+1:dash]
            if model_num in self.best_dict.keys():
                self.best_dict[model_num].append(gen_num)
            else:
                self.best_dict[model_num] = [gen_num]


    def create_flux_maps(self):
        """
        Compute de flux maps. The fluxes are computed on the last
        flux_time_frac portion of the simulations.
        """

        # Copy files to din_dir and add %mod DIN lines.
        for model_num in self.best_dict.keys():
            gen_num = self.best_dict[model_num][0]
            model = "best_{}-{}".format(gen_num, model_num)
            input_path = "{}/{}.ka".format(self.best_dir, model)
            output_path = "{}/{}.ka".format(self.din_dir, model)
            shutil.copyfile(input_path, output_path)
            output_file = open(output_path, "a")
            flux_path = "{}/flux_{}-{}.json".format(self.din_dir,
                                                    gen_num, model_num)
            start_time = self.sim_time * (1-self.flux_time_frac)
            if start_time == 0.0:
                start_time = 0.001
            output_file.write('%mod: alarm {} do $DIN "{}" [true];\n'
                              .format(start_time, flux_path))
            output_file.write('%mod: alarm {} do $DIN "{}" [false];\n'
                              .format(self.sim_time, flux_path))
            output_file.close()
        # Run KaSim with new %mod lines to produce flux maps.
        for model_num in self.best_dict.keys():
            gen_num = self.best_dict[model_num][0]
            model = "best_{}-{}".format(gen_num, model_num)
            input_path = "{}/{}.ka".format(self.din_dir, model)
            output_path = "{}/{}.dat".format(self.din_dir, model)
            command_line = ("KaSim -mode batch --no-log --no-log -u t -p {} "
                            "-l {} -i {} -o {}"
                            .format(self.out_period, self.sim_time,
                                    input_path, output_path))
            os.system(command_line)


    def draw_flux_maps(self):
        """ Draw the flux maps using graphviz. """

        for model_num in self.best_dict.keys():
            # Open flux file.
            #model_num = "867"
            gen_num = self.best_dict[model_num][0]
            flux_path = "{}/flux_{}-{}.json".format(self.din_dir,
                                                    gen_num, model_num)
            flux_file = open(flux_path)
            flux_data = json.load(flux_file)
            din_rules = flux_data["din_rules"]
            din_hits = flux_data["din_hits"]
            din_fluxes = flux_data["din_fluxs"]

            values_path = "{}/values_{}-{}.dat".format(self.din_dir,
                                                    gen_num, model_num)
            values_file = open(values_path, "w")
            values_file.write("                 Source\nTarget          ")
            for i in range(24):
                values_file.write("{:11}|".format(din_rules[i]))
            values_file.write("\n\n")
            for i in range(len(din_fluxes)):
                values_file.write("{:16} ".format(din_rules[i]))
                for j in range(24):
                    value = float(din_fluxes[j][i])
                    if value == 0.0:
                        values_file.write("     0.0    ")
                    else:
                        values_file.write("{:11.6f} ".format(float(din_fluxes[j][i])))
                values_file.write("\n")
          
            # Make a list of base rules.
            base_rules = []
            for rule in din_rules:
                if "Interventions" not in rule and "unbinds" not in rule:
                    space = rule.rfind(" ")
                    base_rule = rule[:space]
                    if "uni" not in base_rule:
                        if base_rule not in base_rules:
                            base_rules.append(base_rule)
            # Get rule indexes.
            rule_indexes = {}
            for base_rule in base_rules:
                rule_indexes[base_rule] = []
                for i in range(len(din_rules)):
                    din_rule = din_rules[i]
                    if "uni" not in din_rule:
                        if base_rule in din_rule:
                            rule_indexes[base_rule].append(i)
            # Get merged hits.
            rule_hits = {}
            hits_list = []
            for base_rule in rule_indexes.keys():
                indexes = rule_indexes[base_rule]
                summ = 0
                for i in indexes:
                    hits = din_hits[i]
                    summ += hits
                #if summ != 0:
                rule_hits[base_rule] = summ
                hits_list.append(summ)
            if len(hits_list) > 0:
                max_hits = max(hits_list)
                min_hits = min(hits_list)
            else:
                max_hits, min_hits = 0, 0
            # Get merged fluxes.
            rule_fluxes = {}
            fluxes_list = []
            for rule1 in base_rules:
                for rule2 in base_rules:
                    edge_name = '"{}" -> "{}"'.format(rule1, rule2)
                    summ = 0
                    for i in rule_indexes[rule1]:
                        for j in rule_indexes[rule2]:
                            flux = din_fluxes[i][j]
                            summ += flux
                    if summ != 0:
                        rule_fluxes[edge_name] = summ
                        fluxes_list.append(abs(summ))
            if len(fluxes_list) > 0:
                max_flux = max(fluxes_list)
            else:
                max_flux = 0

            # Node color range for number of hits.
            num_colors = 7 # Limited by graphviz coloschemes
            palette = "blues8"
            hits_range = max_hits - min_hits
            binwidth = hits_range / num_colors
            # Edge width range for amplitude of flux.
            num_pen = 8
            pen_binwidth = max_flux / num_pen
            # Produce layout.
            radius = 5
            delta = 2*math.pi/len(rule_hits.keys())
            pos_list = []
            for i in range(len(rule_hits.keys())):
                angle = i*delta
                if angle < math.pi:
                    x_sign = 1
                else:
                    y_sign = -1
                if angle < math.pi/2 and angle > 3*math.pi/2:
                    y_sign = 1
                else:
                    y_sign = -1
                x = math.sin(angle)*radius*x_sign
                y = math.cos(angle)*radius*y_sign
                position = "{},{}!".format(x, y)
                pos_list.append(position)

            # Create graph.
            g = graphviz.Digraph(comment="Flux map")
            #g.attr("graph", splines="curved")
            pos_ind = 0
            for rule in rule_hits.keys():
                num_hits = rule_hits[rule]
                if num_hits == min_hits:
                    colorbin = 1
                else:
                    rel_hits = num_hits - min_hits
                    colorbin = math.ceil( rel_hits / binwidth )
                colornum = "{}".format(colorbin)
                g.node(rule, rule, shape="box", style="filled",
                           colorscheme=palette, fillcolor=colornum, 
                           pos=pos_list[pos_ind], fixedsize="true", width="1")
                pos_ind += 1
            for edge in rule_fluxes.keys():
                flux = rule_fluxes[edge]
                flux_str = "{:.4f}".format(flux)
                arrow = edge.index("->")
                source = edge[1:arrow-2]
                target = edge[arrow+4:-1]
                penbin = math.ceil( abs(flux) / pen_binwidth )
                pensize = "{}".format(penbin)
                if flux > 0:
                    edge_color = "green"
                    arrow_shape = "normal"
                else:
                    edge_color = "red"
                    arrow_shape = "tee"
                g.edge(source, target, color=edge_color, arrowhead=arrow_shape,
                       penwidth=pensize, label=flux_str)

            # Write graph to dot file.
            output_path = "{}/flux_{}-{}.dot".format(self.din_dir,
                                                     gen_num, model_num)
            output_file = open(output_path, "w")
            output_file.write(g.source)
            output_file.close()
            # Make a png from the dot file.
            png_path = "{}.png".format(output_path[:-4])
            command_line = ("/usr/bin/neato -Tpng {} > {}".format(output_path,
                                                                png_path))
            os.system(command_line)


    # """""""""""" DIN section End """""""""""""
