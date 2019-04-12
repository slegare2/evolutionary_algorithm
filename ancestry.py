#! /usr/bin/python3

""" Compute ancestral mutation events. """

import os
import shutil


class RuleAncestry:
    """
    Compute ancestral mutation events of a given rule
    from a given Kappa model.
    """

    def __init__(self, kappa_path, rule_id, past_dir, ance_dir):
        """ Initialize RuleAncestry class. """

        self.kappa_path = kappa_path
        self.rule_id = rule_id
        self.past_dir = past_dir
        self.ance_dir = ance_dir
        self.kappa_obj = open(self.kappa_path, "r").readlines()
        # Run class methods.
        self.find_initial_rule()
        self.copy_target_file()
        self.create_outfile()
        self.get_past_models()
        # Main loop.
        self.track_ancestors()


    # ------------- Initialize Section ---------------

    def find_initial_rule(self):
        """ Search for rule_id in the model_file. """

        self.current_model = self.kappa_obj
        self.find_binding_rules()
        rule_found = False
        for i in range(self.bind_start, self.bind_end+1):
            line_str = self.kappa_obj[i][1:]
            quote = line_str.index("'")
            rule_name = line_str[:quote]
            if rule_name == self.rule_id:
                self.rule_line = i
                rule_found = True
                break
        if rule_found == False:
            raise Exception("Rule {} not found in file {}"
                            .format(self.rule_id, self.kappa_path))


    def find_binding_rules(self):
        """ Read input Kappa file and find the range of the binding rules. """

        self.bind_start = 0
        self.bind_end = 0
        for i in range(len(self.current_model)):
            line = self.kappa_obj[i]
            if "// Binary binding rules." in line:
                self.bind_start = i + 1
            if "// Unary binding rules." in line:
                self.bind_end = i - 1


    def copy_target_file(self):
        """ Copy the model_file to ance_dir. """

        if "/" in self.kappa_path:
            slash = self.kappa_path.rfind("/")
            self.kappa_file = self.kappa_path[slash+1:]
        else:
            self.kappa_file = self.kappa_path
        to_path = "{}/{}".format(self.ance_dir, self.kappa_file)
        shutil.copyfile(self.kappa_path, to_path)


    def create_outfile(self):
        """ Create file where to write ancestry results. """ 

        output_path = "{}/{}.txt".format(self.ance_dir, self.kappa_file[:-3])
        self.output_file = open(output_path, "w")


    def get_past_models(self):
        """ Get the list of all past generation models. """

        file_list = os.listdir(self.past_dir)
        self.past_list = []
        for f in file_list:
            dash = f.index("-")
            generation = int(f[4:dash])
            model_id = int(f[dash+1:-3])
            self.past_list.append([generation, model_id])

    # ------------- Initialize Section End ---------------


    def track_ancestors(self):
        """ Track mutations in the past of given rule from given model. """

        self.root_reached = False
        self.previous_model = self.kappa_obj
        self.current_model = self.kappa_obj
        while self.root_reached == False:
            self.find_target_rule()
            if self.mutation_found == True:

            self.previous_model = self.current_model 
            self.fetch_father()
            if self.father_id == 0:
                self.root_reached = True


    def find_target_rule(self):
        """
        Find target rule in model file. If the rule is not found, it has been
        mutated and its precedent form should be found instead.
        """

        self.find_binding_rules()
        self.rule_found = False
        self.mutation_found = False
        for i in range(self.bind_start, self.bind_end+1):
            line_str = self.current_model[i][1:]
            quote = line_str.index("'")
            rule_name = line_str[:quote]
            if rule_name == self.rule_id:
                self.rule_line = i
                self.rule_found = True
                break
        if rule_found == False:
            self.rule_id = self.rule_id[:-1]
        for i in range(self.bind_start, self.bind_end+1):
            line_str = self.current_model[i][1:]
            quote = line_str.index("'")
            rule_name = line_str[:quote]
            if rule_name == self.rule_id:
                self.rule_line = i
                self.mutation_found = True
                break


    def fetch_father(self):
        """
        Get the father model from directory past_dir. If the model was
        preserved for many generations, take the file from the first
        generation where it appeared.
        """

        self.get_father_id()
        father_list = []
        for past_model in self.past_list:
            if past_model[1] == self.father_id:
                father_list.append(past_model)
        sorted_list = sorted(father_list, key=lambda x: x[0])
        self.father_file = ("gen_{}-{}.ka".format(sorted_list[0][0],
                                                  sorted_list[0][1]))
        self.father_path = "{}/{}".format(self.past_dir, self.father_file)
        self.current_model = open(self.father_path, "r").readlines()


    def get_father_id(self):
        """ Read the father id from the model. """

        tokens = self.kappa_obj[0].split()
        father = tokens[-1]
        if "-" in father:
            dash = father.index("-")
            self.father_id = int(father[dash+1:-3])
        else:
            self.father_id = 0


class ModelAncestry:
    """ Compute ancestral mutation events of a given Kappa model. """

    def __init__(self, model_file, past_dir):
        """ Initialize RuleAncestry class. """

        self.model_file = model_file
        self.past_dir = past_dir

