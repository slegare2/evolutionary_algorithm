#! /usr/bin/python3

""" Compute ancestral mutation events. """

import os
import math
import shutil
import graphviz


class Ancestry:
    """ Compute all the ancestry of an evolutionary algorithm run. """

    def __init__(self, fit_dir, past_dir, ance_dir):
        """ Initialize Ancestry class. """

        self.fit_dir = fit_dir
        self.past_dir = past_dir
        self.ance_dir = ance_dir
        self.nodes = []
        self.edges = []
        self.min_fitness = 100000
        self.max_fitness = -100000
        self.max_diff = -100000
        self.mutations = {}
        self.fitness_diffs = {}
        self.fathers = {}
        self.sons = {}
        self.num_rules = {}
        # Run class methods.
        self.compute_ancestry()
        self.draw_ancestry()


    def compute_ancestry(self):
        """ Main method. """

        fitness_files = self.get_files(self.fit_dir)
        self.last_gen = len(fitness_files)
        self.get_node_zero()
        for gen in range(1, self.last_gen+1):
            self.models_fitness = {}
            self.add_nodes(gen)
            self.add_survivor_edges(gen)
            self.add_mutation_edges(gen)
        best_num = self.find_best()
        self.list_parents([best_num])
        self.center_alpha_models()


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
            number = int(f[dash+1:dot])
            file_dicts.append({"file": name, "num": number})
        sorted_dicts = sorted(file_dicts, key=lambda x: x["num"])
        file_list = []
        for d in sorted_dicts:
            file_list.append(d["file"])
        return file_list


    def get_node_zero(self):
        """ Find the fitness of model 0."""

        input_path = "{}/fitness_gen-1.txt".format(self.fit_dir)
        input_file = open(input_path).readlines()
        for i in range(2, len(input_file)):
            tokens = input_file[i].split()
            model = tokens[0]
            dash = model.index("-")
            model_num = model[dash+1:]
            if model_num == "0":
                fitness_zero = float(tokens[-1])
                nrules_zero = tokens[1]
                self.nodes.append([[model, nrules_zero, fitness_zero]])
                break


    def add_nodes(self, generation):
        """ Get the nodes (the models) of a given generation. """

        input_path = "{}/fitness_gen-{}.txt".format(self.fit_dir, generation)
        input_file = open(input_path).readlines()
        self.current_models = []
        rank_nodes = []
        for i in range(2, len(input_file)):
            tokens = input_file[i].split()
            model = tokens[0]
            nrules = tokens[1]
            fitness = float(tokens[-1])
            if fitness < self.min_fitness:
                self.min_fitness = fitness
            if fitness > self.max_fitness:
                self.max_fitness = fitness
            rank_nodes.append([model, nrules, fitness])
            self.current_models.append(model)
            self.models_fitness[model] = fitness
        self.nodes.append(rank_nodes)


    def add_survivor_edges(self, generation):
        """ Get the edges (the mutations) of a given generation. """

        if generation == 1:
            input_file = ["", "", "abc-0 0.0"]
        else:
            input_path = ("{}/fitness_gen-{}.txt"
                          .format(self.fit_dir, generation-1))
            input_file = open(input_path).readlines()
        self.survivors = []
        self.rank_edges = []
        for i in range(2, len(input_file)):
            tokens = input_file[i].split()
            prev_model = tokens[0]
            if prev_model in self.current_models:
                self.rank_edges.append([prev_model, prev_model, 0])
                self.survivors.append(prev_model)
                dash = prev_model.rfind("-")
                model_num = prev_model[dash+1:]
                node_id = "node_{}_{}".format(generation, model_num)
                self.mutations[node_id] = "Clone"
                prev_id = "node_{}_{}".format(generation-1, model_num)
                if prev_id not in self.sons.keys():
                    self.sons[prev_id] = [prev_model]
                else:
                    self.sons[prev_id].append(prev_model)
            else:
                fitness = float(tokens[-1])
                self.models_fitness[prev_model] = fitness


    def add_mutation_edges(self, generation):
        """ Get the edges (the mutations) of a given generation. """

        input_path = "{}/fitness_gen-{}.txt".format(self.fit_dir, generation)
        input_file = open(input_path).readlines()
        for i in range(2, len(input_file)):
            tokens = input_file[i].split()
            model = tokens[0]
            if model not in self.survivors:
                dash = model.rfind("-")
                model_num = int(model[dash+1:])
                model_path = ("{}/gen_{}-{}.ka"
                              .format(self.past_dir, generation, model_num))
                model_file = open(model_path).readline()
                father = model_file.split()[2][:-3]
                curr_fitness = self.models_fitness[model]
                prev_fitness = self.models_fitness[father]
                fitness_diff = curr_fitness - prev_fitness
                if fitness_diff > self.max_diff:
                    self.max_diff = fitness_diff
                curr_id = "node_{}_{}".format(generation, model_num)
                self.rank_edges.append([father, model, fitness_diff])
                self.retrace_mutation(father, model_path, curr_id, generation)
                self.fitness_diffs[curr_id] = fitness_diff
                self.fathers[curr_id] = [father, prev_fitness]
                father_dash = father.rfind("-")
                father_num = int(father[dash+1:])
                prev_id = "node_{}_{}".format(generation-1, father_num)
                if prev_id not in self.sons.keys():
                    self.sons[prev_id] = [model]
                else:
                    self.sons[prev_id].append(model)
        self.edges.append(self.rank_edges)


    def retrace_mutation(self, father, son_path, son_id, generation):
        """ Find what mutation occured between a father a son model. """

        mutation_str = ""
        dash = father.rfind("-")
        father_num = int(father[dash+1:])
        if generation == 1:
            father_path = ("{}/gen_{}-{}.ka"
                           .format(self.past_dir, generation, father_num))
        else:
            father_path = ("{}/gen_{}-{}.ka"
                           .format(self.past_dir, generation-1, father_num))
        father_file = open(father_path).readlines()
        son_file = open(son_path).readlines()
        father_start, father_end = self.find_binding_rules(father_file)
        son_start, son_end = self.find_binding_rules(son_file)
        n_lines_father = father_end - father_start
        added_lines = 0
        for i in range(n_lines_father+1):
            father_line = father_file[i + father_start]
            son_line = son_file[i + son_start + added_lines]
            if "// Unary binding rules." not in father_line:
                if father_line != son_line:
                    # Check if mutation is a refinement or a change of rate.
                    father_quote = father_line[1:].index("'")+1
                    father_rname = father_line[1:father_quote]
                    son_quote = son_line[1:].index("'")+1
                    son_rname = son_line[1:son_quote]
                    if son_rname != father_rname:
                        mutation_type = "refinement"
                        # Find how many rules were added from the refinement.
                        n_refine = 0
                        for j in range(i+son_start+added_lines,
                                       father_end+added_lines):
                            check_line = son_file[j]
                            if "// Unary binding rules." in check_line:
                                break
                            check_quote = check_line[1:].index("'")+1
                            check_rname = check_line[1:check_quote]
                            if father_rname in check_rname:
                                n_refine += 1
                            else:
                                break
                        if "uni" not in father_rname:
                            mutation_str += father_line
                            mutation_str += "        ||\n        V\n"
                            for j in range(n_refine):
                                new_line = son_file[i + son_start + added_lines + j]
                                mutation_str += new_line
                            mutation_str += "\n"
                        added_lines = added_lines + n_refine - 1
                    else:
                       mutation_type = "rate change"
                       mutation_str += father_line
                       mutation_str += "        ||\n        V\n"
                       mutation_str += son_line
                       mutation_str += "\n"
        self.mutations[son_id] = [mutation_type, mutation_str]


    def find_binding_rules(self, model_file):
        """ Read input Kappa file and find the range of the binding rules. """

        bind_start = 0
        bind_end = 0
        for i in range(len(model_file)):
            line = model_file[i]
            if "// Binary binding rules." in line:
                bind_start = i + 1
            if "// Unbinding rules." in line:
                bind_end = i - 1

        return bind_start, bind_end


    def find_best(self):
        """ Find the model with highest fitness in the last generation. """

        input_path = "{}/fitness_gen-{}.txt".format(self.fit_dir,
                                                    self.last_gen)
        input_file = open(input_path, "r").readlines()
        tokens = input_file[2].split()
        best_model = tokens[0]
        dash = best_model.index("-")
        best_num = int(best_model[dash+1:])

        return best_num


    def list_parents(self, model_nums):
        """ Make a list of the ancestors of a given list of model numbers. """

        past_files = self.get_files(self.past_dir)
        past_list = []
        for f in past_files:
            dash = f.rfind("-")
            underscore = f.rfind("_")
            f_gen = int(f[underscore+1:dash])
            f_num = int(f[dash+1:])
            past_list.append([f_gen, f_num])
        self.parents_highlight = ["node_0_0"]
        for model_num in model_nums:
            self_num = model_num
            for gen in range(self.last_gen, 0, -1):
                # Search for the model itself in current generation.
                model_found = False
                for past_mod in past_list:
                    if past_mod[0] == gen and past_mod[1] == self_num:
                        node_id = "node_{}_{}".format(gen, self_num)
                        self.parents_highlight.append(node_id)
                        model_found = True
                        break
                # If the model itself is not found, search for its father.
                if model_found == False:
                    # First, find the model itself in subsequent generation.
                    for past_mod in past_list:
                        if past_mod[0] == gen+1 and past_mod[1] == self_num:
                            model_path = ("{}/gen_{}-{}.ka".format(
                                          self.past_dir, gen+1, self_num))
                            model_file = open(model_path, "r").readline()
                            father = model_file.split()[2][:-3]
                            father_dash = father.rfind("-")
                            father_num = int(father[father_dash+1:])
                            break
                    # Then, search for the father in current generation.
                    for past_mod in past_list:
                        if past_mod[0] == gen and past_mod[1] == father_num:
                            node_id = "node_{}_{}".format(gen, father_num)
                            self.parents_highlight.append(node_id)
                            self_num = father_num
                            break


    def center_alpha_models(self):
        """ Reorder nodes to center models with most sons. """

        self.centered_nodes = []
        for gen in range(0, self.last_gen+1):
            rank_nodes = self.nodes[gen]
            tmp_list = []
            for node in rank_nodes:
                model = node[0]
                dash = model.index("-")
                model_num = model[dash+1:]
                node_id = "node_{}_{}".format(gen, model_num)
                try:
                    sons_list = self.sons[node_id]
                    n_sons = len(sons_list)
                except:
                    n_sons = 0
                tmp_node = node.copy()
                tmp_node.append(n_sons)
                tmp_list.append(tmp_node)
            sorted_list = sorted(tmp_list, key=lambda x: x[3], reverse=True)
            centered_list = []
            for i in range(0, len(sorted_list)):
                if gen%2 == 0:
                    if i%2 != 0:
                        centered_list.append(sorted_list[i])
                    else:
                        centered_list.insert(0, sorted_list[i])
                else:
                    if i%2 != 0:
                        centered_list.insert(0, sorted_list[i])
                    else:
                        centered_list.append(sorted_list[i])
            self.centered_nodes.append(centered_list)



    def draw_ancestry(self):
        """ Output the ancestry in dot format using graphviz. """

        # Custom position of node if using neato (ignored if using dot).
        x_spacing = 1.5
        y_spacing = 3
        y = len(self.centered_nodes) * y_spacing
        # Node color range for fitness.
        num_colors = 8 # Limited by graphviz coloschemes
        colorschemes = "blues"
        palette = "{}{}".format(colorschemes, num_colors)
        fitness_range = self.max_fitness - self.min_fitness
        binwidth = fitness_range / num_colors
        # Edge width range for fitness difference between models.
        num_pen = 8
        pen_binwidth = self.max_diff / num_pen
        # Build graph.
        g = graphviz.Graph(comment="Model ancestry")
        for gen in range(0, len(self.centered_nodes)):
        #for gen in range(0, 150):
            x = 0
            generation_nodes = self.centered_nodes[gen]
            for node in generation_nodes:
                dash = node[0].index("-")
                node_num = node[0][dash+1:]
                node_id = "node_{}_{}".format(gen, node_num)
                nrules = node[1]
                fitness = node[2]
                if fitness == self.min_fitness:
                    colorbin = 1
                else:
                    rel_fit = fitness - self.min_fitness
                    colorbin = math.ceil( rel_fit / binwidth )
                colornum = "{}".format(colorbin)
                position = "{},{}!".format(x, y)
                if gen > 0:
                    if self.mutations[node_id] != "Clone":
                        clickmsg =  ("\nThis is model {} at generation {}\n"
                                     .format(node[0], gen))
                        clickmsg += ("Number of rules = {}\n".format(nrules))
                        fitness_diff = self.fitness_diffs[node_id]
                        father = self.fathers[node_id][0]
                        if fitness_diff >= 0:
                            sign = "+"
                        else:
                            sign = ""
                        clickmsg +=  ("Fitness = {:.3f}  ({}{:.3f})\n\n"
                                      .format(fitness, sign, fitness_diff))
                        clickmsg +=  ("Mutation compared to previous "
                                      "model {} :\n\n".format(father))
                        #clickmsg += ("({})\n"
                        #             .format(self.mutations[node_id][0]))
                        clickmsg += ("{}".format(self.mutations[node_id][1]))
                    else:
                        clickmsg = ("\nClone of model {} at generation {}.\n\n"
                                    .format(node[0], gen))
                        clickmsg += "Fitness = {}".format(fitness)
                else:
                    clickmsg = ("\nStarting model {} ."
                                .format(node[0]))
                if node_id in self.parents_highlight:
                    pensize = "2"
                    bordercolor = "#CD5C5C"
                    nodeshape = "doublecircle"
                else:
                    pensize = "1"
                    bordercolor = "#000000"
                    nodeshape = "circle"
                g.node(node_id, node_num, shape=nodeshape, style="filled",
                       colorscheme=palette, fillcolor=colornum, pos=position,
                       URL=clickmsg, penwidth=pensize, fixedsize="true",
                       width="0.6", color=bordercolor)
                x += x_spacing
            y += -y_spacing
        for gen in range(1, len(self.edges)+1):
        #for gen in range(1, 150):
            generation_edges = self.edges[gen-1]
            for edge in generation_edges:
                dash1 = edge[0].index("-")
                node1_num = edge[0][dash+1:]
                node1_id = "node_{}_{}".format(gen-1, node1_num)
                dash2 = edge[1].index("-")
                node2_num = edge[1][dash+1:]
                node2_id = "node_{}_{}".format(gen, node2_num)
                fit_diff = edge[2]
                if node1_num == node2_num:
                    g.attr("edge", style="dashed", color="grey", penwidth="3")
                else:
                    if fit_diff <= 0:
                        penbin = 1
                    else:
                        penbin = math.ceil( fit_diff / pen_binwidth )
                    pensize = "{}".format(penbin)
                    if self.mutations[node2_id][0] == "refinement":
                        edge_color = "green"
                    else:
                        edge_color = "black"
                    g.attr("edge", style="solid", color=edge_color, penwidth=pensize)
                g.edge(node1_id, node2_id)
        outfile = open("ancestry.dot", "w")
        outfile.write(g.source)
        outfile.close()
        

#class RuleAncestry:
#    """
#    Compute ancestral mutation events of a given rule
#    from a given Kappa model.
#    """
#
#    def __init__(self, kappa_path, rule_id, past_dir, ance_dir):
#        """ Initialize RuleAncestry class. """
#
#        self.kappa_path = kappa_path
#        self.rule_id = rule_id
#        self.past_dir = past_dir
#        self.ance_dir = ance_dir
#        self.kappa_model = open(self.kappa_path, "r").readlines()
#        # Run class methods.
#        self.check_rule()
#        self.copy_target_file()
#        self.create_outfile()
#        self.get_past_models()
#        # Main loop.
#        self.track_ancestors()
#
#
#    # ------------- Initialize Section ---------------
#
#    def check_rule(self):
#        """ Search for rule_id in the model_file. """
#
#        self.find_binding_rules(self.kappa_model)
#        rule_found = False
#        for i in range(self.bind_start, self.bind_end+1):
#            line_str = self.kappa_model[i][1:]
#            quote = line_str.index("'")
#            rule_name = line_str[:quote]
#            if rule_name == self.rule_id:
#                self.rule_line = i
#                rule_found = True
#                break
#        if rule_found == False:
#            raise Exception("Rule {} not found in file {}"
#                            .format(self.rule_id, self.kappa_path))
#
#
#    def find_binding_rules(self, model_file):
#        """ Read input Kappa file and find the range of the binding rules. """
#
#        self.bind_start = 0
#        self.bind_end = 0
#        for i in range(len(model_file)):
#            line = model_file[i]
#            if "// Binary binding rules." in line:
#                self.bind_start = i + 1
#            if "// Unary binding rules." in line:
#                self.bind_end = i - 1
#
#
#    def copy_target_file(self):
#        """ Copy the model_file to ance_dir. """
#
#        if "/" in self.kappa_path:
#            slash = self.kappa_path.rfind("/")
#            self.kappa_file = self.kappa_path[slash+1:]
#        else:
#            self.kappa_file = self.kappa_path
#        to_path = "{}/{}".format(self.ance_dir, self.kappa_file)
#        shutil.copyfile(self.kappa_path, to_path)
#
#
#    def create_outfile(self):
#        """ Create file where to write ancestry results. """ 
#
#        output_path = "{}/{}.txt".format(self.ance_dir, self.kappa_file[:-3])
#        self.output_file = open(output_path, "w")
#
#
#    def get_past_models(self):
#        """ Get the list of all past generation models. """
#
#        file_list = os.listdir(self.past_dir)
#        self.past_list = []
#        for f in file_list:
#            dash = f.index("-")
#            generation = int(f[4:dash])
#            model_id = int(f[dash+1:-3])
#            self.past_list.append([generation, model_id])
#
#    # ------------- Initialize Section End ---------------
#
#
#    def track_ancestors(self):
#        """ Track mutations in the past of given rule from given model. """
#
#        root_reached = False
#        current_file = self.kappa_path
#        current_rule = self.rule_id
#        self.out_str = "{}\n".format(self.kappa_file)
#        while root_reached == False:
#            current_model = open(current_file, "r").readlines()
#            current_rule, current_line = self.find_target_rule(current_rule,
#                                                               current_model)
#            next_file, father_id = self.find_father(current_model)
#            next_model = open(next_file, "r").readlines()
#            next_rule, next_line = self.find_target_rule(current_rule,
#                                                         next_model)
#            if next_rule != current_rule: # A mutation was found.        
#                self.out_str += "|\n"
#                slash = current_file.rfind("/")
#                self.out_str += "{:20} ".format(current_file[slash+1:])
#                self.out_str += "{}".format(current_model[current_line])
#                slash = next_file.rfind("/")
#                self.out_str += "{:20} ".format(next_file[slash+1:])
#                self.out_str += "{}".format(next_model[next_line])
#            else:
#                pass
#            current_file = next_file
#            if father_id == 0:
#                root_reached = True
#        self.output_file.write("{}".format(self.out_str))
#
#
#    def find_father(self, model_file):
#        """
#        Get the father model from directory past_dir. If the model was
#        preserved for many generations, take the file from the first
#        generation where it appeared.
#        """
#
#        tokens = model_file[0].split()
#        father = tokens[-1]
#        if "-" in father:
#            dash = father.index("-")
#            father_id = int(father[dash+1:-3])
#        else:
#            father_id = 0
#        father_list = []
#        for past_model in self.past_list:
#            if past_model[1] == father_id:
#                father_list.append(past_model)
#        sorted_list = sorted(father_list, key=lambda x: x[0])
#        father_gen = sorted_list[0][0]
#        father_file = "gen_{}-{}.ka".format(father_gen, father_id)
#        father_path = "{}/{}".format(self.past_dir, father_file)
#
#        return father_path, father_id
#
#
#    def find_target_rule(self, rule, model_file):
#        """
#        Find target rule in model file. If the rule is not found, it has been
#        mutated and its precedent form should be found instead.
#        """
#
#        self.find_binding_rules(model_file)
#        rule_found = False
#        for i in range(self.bind_start, self.bind_end+1):
#            line_str = model_file[i][1:]
#            quote = line_str.index("'")
#            rule_name = line_str[:quote]
#            if rule_name == rule:
#                rule_line = i
#                rule_found = True
#                break
#        if rule_found == False:
#            prev_rule = rule[:-1]
#            mutation_found = True
#            for i in range(self.bind_start, self.bind_end+1):
#                line_str = model_file[i][1:]
#                quote = line_str.index("'")
#                rule_name = line_str[:quote]
#                if rule_name == prev_rule:
#                    rule_line = i
#                    mutation_found = True
#                    break
#            if mutation_found == False:
#                raise Exception("Could not find rule {} or {}."
#                               .format(rule, prev_rule))
#            else:
#                rule = prev_rule
#
#        return rule, rule_line
#
#
#class ModelAncestry:
#    """ Compute ancestral mutation events of a given Kappa model. """
#
#    def __init__(self, model_file, past_dir):
#        """ Initialize RuleAncestry class. """
#
#        self.model_file = model_file
#        self.past_dir = past_dir

