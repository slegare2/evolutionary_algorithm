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
        self.kappa_model = open(self.kappa_path, "r").readlines()
        # Run class methods.
        self.check_rule()
        self.copy_target_file()
        self.create_outfile()
        self.get_past_models()
        # Main loop.
        self.track_ancestors()


    # ------------- Initialize Section ---------------

    def check_rule(self):
        """ Search for rule_id in the model_file. """

        self.find_binding_rules(self.kappa_model)
        rule_found = False
        for i in range(self.bind_start, self.bind_end+1):
            line_str = self.kappa_model[i][1:]
            quote = line_str.index("'")
            rule_name = line_str[:quote]
            if rule_name == self.rule_id:
                self.rule_line = i
                rule_found = True
                break
        if rule_found == False:
            raise Exception("Rule {} not found in file {}"
                            .format(self.rule_id, self.kappa_path))


    def find_binding_rules(self, model_file):
        """ Read input Kappa file and find the range of the binding rules. """

        self.bind_start = 0
        self.bind_end = 0
        for i in range(len(model_file)):
            line = model_file[i]
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

        root_reached = False
        current_file = self.kappa_path
        current_rule = self.rule_id
        self.out_str = "{}\n".format(self.kappa_file)
        while root_reached == False:
            current_model = open(current_file, "r").readlines()
            current_rule, current_line = self.find_target_rule(current_rule,
                                                               current_model)
            next_file, father_id = self.find_father(current_model)
            next_model = open(next_file, "r").readlines()
            next_rule, next_line = self.find_target_rule(current_rule,
                                                         next_model)
            if next_rule != current_rule: # A mutation was found.        
                self.out_str += "|\n"
                slash = current_file.rfind("/")
                self.out_str += "{:20} ".format(current_file[slash+1:])
                self.out_str += "{}".format(current_model[current_line])
                slash = next_file.rfind("/")
                self.out_str += "{:20} ".format(next_file[slash+1:])
                self.out_str += "{}".format(next_model[next_line])
            else:
                pass
            current_file = next_file
            if father_id == 0:
                root_reached = True
        self.output_file.write("{}".format(self.out_str))


    def find_father(self, model_file):
        """
        Get the father model from directory past_dir. If the model was
        preserved for many generations, take the file from the first
        generation where it appeared.
        """

        tokens = model_file[0].split()
        father = tokens[-1]
        if "-" in father:
            dash = father.index("-")
            father_id = int(father[dash+1:-3])
        else:
            father_id = 0
        father_list = []
        for past_model in self.past_list:
            if past_model[1] == father_id:
                father_list.append(past_model)
        sorted_list = sorted(father_list, key=lambda x: x[0])
        father_gen = sorted_list[0][0]
        father_file = "gen_{}-{}.ka".format(father_gen, father_id)
        father_path = "{}/{}".format(self.past_dir, father_file)

        return father_path, father_id


    def find_target_rule(self, rule, model_file):
        """
        Find target rule in model file. If the rule is not found, it has been
        mutated and its precedent form should be found instead.
        """

        self.find_binding_rules(model_file)
        rule_found = False
        for i in range(self.bind_start, self.bind_end+1):
            line_str = model_file[i][1:]
            quote = line_str.index("'")
            rule_name = line_str[:quote]
            if rule_name == rule:
                rule_line = i
                rule_found = True
                break
        if rule_found == False:
            prev_rule = rule[:-1]
            mutation_found = True
            for i in range(self.bind_start, self.bind_end+1):
                line_str = model_file[i][1:]
                quote = line_str.index("'")
                rule_name = line_str[:quote]
                if rule_name == prev_rule:
                    rule_line = i
                    mutation_found = True
                    break
            if mutation_found == False:
                raise Exception("Could not find rule {} or {}."
                               .format(rule, prev_rule))
            else:
                rule = prev_rule

        return rule, rule_line


class ModelAncestry:
    """ Compute ancestral mutation events of a given Kappa model. """

    def __init__(self, model_file, past_dir):
        """ Initialize RuleAncestry class. """

        self.model_file = model_file
        self.past_dir = past_dir

