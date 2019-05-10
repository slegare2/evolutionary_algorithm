#! /usr/bin/python3

"""
Selects a random binding rule from a kappa file and refine it by adding one
condition from the possible bindings found in a base interaction matrix.

Example usage:

# User parameter. Maximum number of edges from the binding partner.
max_bonds = 3

# Input Kappa file and interaction matrix.
kappa_file = "random_inters.ka"
matrix_file = "inter_matrix.txt"

# Use class MutateRule to mutate a random rule.
mutated_model = MutateRule(kappa_file, matrix_file, max_bonds)
print(str(mutated_model))

# Write the mutated model to output file.
out_file = outfilename(kappa_file)
out_obj = open(out_file, "w")
out_obj.write(repr(mutated_model))
out_obj.close()

"""

import random
import string


class MutateRule:
    """
    Mutate a random rule from a Kappa model into a set of refined rules.
    Possible bonds to add in refinements are taken from the interaction matrix
    file. New bonds added can form rings if the added partner type corresponds
    to an agent that is already present in the initial rule.

    To test a mutation on a specfic rule, put the desired rule into a
    kappa_file and provide select_line, the line of the rule in the file.
    A specific mutation_site can also be given as A(x). To refer to, let's
    say, the second agent A of a rule, write A-2(x).
    """

    def __init__(self, kappa_file, matrix_file, max_bnd,
                 select_line=None, mutation_site=None):
        """ Initialize MutateRule class. """

        self.kappa_file = kappa_file
        self.matrix_file = matrix_file
        self.max_bnd = max_bnd
        self.select_line = select_line
        self.mutation_site = mutation_site
        self.kappa_obj = open(self.kappa_file, "r").readlines()
        self.matrix_obj = open(self.matrix_file, "r").readlines()
        # Run class methods.
        self.random_select()
        self.build_new_rules()
        self.build_new_rings()
        self.build_new_file()


    def read_inter_matrix(self):
        """
        Read pairwise interaction matrix from file into a list of lists.
        """

        self.inter_matrix = []
        self.protein_ids = self.matrix_obj[6].split()
        for i in range(8, len(self.matrix_obj)):
            tokens = self.matrix_obj[i].split()
            vector = tokens[1:]
            row = []
            for str_value in vector:
                row.append(int(str_value))
            self.inter_matrix.append(row)
        #print("Interaction matrix:")
        #for row in inter_matrix:
        #    print(row)


    def find_binding_rules(self):
        """ Read input Kappa file and find the range of the binding rules. """
    
        self.bind_start = 0
        self.bind_end = 0
        for i in range(len(self.kappa_obj)):
            line = self.kappa_obj[i]
            if "// Binary binding rules." in line:
                self.bind_start = i + 1
            if "// Unary binding rules." in line:
                self.bind_end = i - 1


    def select_rule(self):
        """
        Select a binding rule, either randomly or from user specified
        select_line. A binding rule should be present at the specified line.
        """

        self.select_mechanism()
        if self.select_line == None: # Select randomly.
            possible_lines = []
            for rule_obj in self.rule_lines:
                if rule_obj["mechanism"] == self.selected_mechanism:
                    possible_lines.append(rule_obj["line"])
            self.random_line = random.choice(possible_lines)
            self.selected_rule = self.kappa_obj[self.random_line][:-1]
        else: # Get the rule from specified line.
            self.selected_rule = self.kappa_obj[self.select_line][:-1]
            self.random_line = self.select_line


    def select_mechanism(self):
        """
        Select a mechanism, either randomly or from from user specified
        select_line. A mechanism is for example 'A binds B', with all its
        possible refinements.
        """

        self.rule_lines = []
        mechanism_list = []
        for i in range(self.bind_start, self.bind_end+1):
            rule_str = self.kappa_obj[i][:-1]
            quote = rule_str.rfind("'")
            rule_name = rule_str[1:quote]
            space = rule_name.rfind(" ")
            mechanism = rule_name[:space]
            if "cls" in mechanism:
                pos = mechanism.index("cls")
                mechanism = mechanism[:pos-1]
            rule_obj = {"mechanism": mechanism,
                        "rule": rule_str,
                        "line": i}
            self.rule_lines.append(rule_obj)
            if mechanism not in mechanism_list:
                mechanism_list.append(mechanism)
            if i == self.select_line:
                self.selected_mechanism = mechanism
        if self.select_line == None:
            self.selected_mechanism = random.choice(mechanism_list)


    def rule_name_rate(self):
        """ Find the name and rate of the selected rule. """ 

        self.name_end = self.selected_rule.rfind("'")
        self.rule_name = self.selected_rule[1:self.name_end]
        self.rate_begin = self.selected_rule.rfind("@")
        if "//" in self.selected_rule:
            slashes = self.selected_rule.index("//")
            self.rule_rate = self.selected_rule[self.rate_begin+2:slashes-1]
        else:
            self.rule_rate = self.selected_rule[self.rate_begin+2:]
        # Get info on the corresponding unary rule.
        self.find_unary()
        self.unary_begin = self.unary_rule.rfind("@")
        if "//" in self.unary_rule:
            slashes = self.unary_rule.index("//")
            self.unary_rate = self.unary_rule[self.unary_begin+2:slashes-1]
        else:
            self.unary_rate = self.unary_rule[self.unary_begin+2:]


    def find_unary(self):
        """
        Find the unimolecular rule that corresponds to the selected rule.
        """
        
        tokens = self.rule_name.split()
        try:
            rule_numsuf = int(tokens[-1])
            space = self.rule_name.rfind(" ")
            rule_pref = self.rule_name[:space]
        except:
            rule_numsuf = ""
        self.unary_name = "{} uni {}".format(self.selected_mechanism,
                                             rule_numsuf)
        for i in range(self.bind_end+1, len(self.kappa_obj)):
            line = self.kappa_obj[i]
            if self.unary_name in line:
                if line [:2] == "//":
                   self.unary_commented = True
                   self.unary_rule = line[2:-1]
                else:
                   self.unary_commented = False
                   self.unary_rule = line[:-1]
                self.unary_line = i
                break


    def agents_connections(self):
        """
        Builds two dictionaries and two lists.
        agents_dictionary contains the agents and sites of the rule.
        connections_dictionary contains only the bonds (excluding "/"s).
        If an agent of the same type is found more than once in the rule, a
        numerical suffix is added to distinguish them in the dictionaries.
        The first list contains the agents of the rule in the form of a list
        rather than a dictionary. The second list contains the two partners
        to be bound by the rule.
        """
    
        self.agents_dictionary = {}
        self.connections_dictionary = {}
        self.partners_list = []
        agents_str = self.selected_rule[self.name_end+2:self.rate_begin-1]
        self.agents_list = agents_str.split(', ')
        seen_ids = {}
        for agent in self.agents_list:
            parenthesis = agent.index("(")
            ptn_type = agent[:parenthesis]
            ptn_id = ptn_type
            if ptn_id in seen_ids.keys():
                number = seen_ids[ptn_id] + 1
                seen_ids[ptn_id] = seen_ids[ptn_id] + 1
                ptn_id = "{}-{}".format(ptn_id, number)
            else:
                seen_ids[ptn_id] = 1
            sites = agent[parenthesis+1:-1].split()
            site_dict = {}
            conn_dict = {}
            for site in sites:
                bracket = site.index("[")
                site_id = site[:bracket]
                bnd_state = site[bracket+1:-1]
                site_dict[site_id] = bnd_state
                if "/" in bnd_state:
                    self.partners_list.append(ptn_id)
                if "/" not in bnd_state and "." not in bnd_state:
                    conn_dict[site_id] = bnd_state
            self.agents_dictionary[ptn_id] = {"type": ptn_type,
                                              "loc": seen_ids[ptn_type],
                                              "sites": site_dict}
            if len(conn_dict) != 0:
                self.connections_dictionary[ptn_id] = {"type": ptn_type,
                                                       "loc": seen_ids[ptn_type],
                                                       "connect": conn_dict}
        if len(self.partners_list) != 2:
            raise Exception("There should be only one bond formed"
                            "by any rule in the current context.")


    def select_partner(self):
        """
        Randomly select one of the binding partners to hold the refinement.
        """
    
        self.random_partner = random.randrange(2)
        self.selected_partner = self.partners_list[self.random_partner]
        other_partner = 0
        if self.random_partner == 0:
            other_partner = 1
        self.unselected_partner = self.partners_list[other_partner]


    def find_target(self, src_agent, src_site, src_bnd):
        """ Find the binding partner of a given agent site. """
    
        bond_found = False
        targ_agents = []
        for agent in self.connections_dictionary.keys():
            if agent != src_agent:
                potential_targs = self.connections_dictionary[agent]["connect"]
                for site in potential_targs.keys():
                    bnd_num = potential_targs[site]
                    if bnd_num == src_bnd:
                        targ_agents.append(agent)
                        #print("Bond {}({}[{}]), {}({}[{}]) found"
                        #      .format(src_agent, src_site,
                        #              src_bnd, agent,
                        #              site, bnd_num))
                        bond_found = True
                        break
        if bond_found == False:
            raise Exception("Target of {}({}[{}]) not found."
                            .format(src_agent, src_site,
                                    bnd_num))    
        return targ_agents


    def find_complex(self, sel_partner):
        """ 
        Find all the agents already connected to the selected_partner
        using connections_dict.
        """
    
        connected_agents = [sel_partner]
        source_agents = [sel_partner]
        seen_bnd_ids = []
        num_bonds = 1
        while len(source_agents) > 0 and num_bonds <= self.max_bnd:
            target_agents = []
            num_bonds += 1
            for source_agent in source_agents:
                try:
                    source_sites = (self.
                        connections_dictionary[source_agent]["connect"])
                    for source_site in source_sites.keys():
                        source_bnd = source_sites[source_site]
                        if source_bnd not in seen_bnd_ids:
                            seen_bnd_ids.append(source_bnd)
                            target_agents += self.find_target(source_agent,
                                source_site, source_bnd)
                except:
                    pass
            source_agents = target_agents
            connected_agents += target_agents
        connected_set = []
        for agent in connected_agents:
            if agent not in connected_set:
                connected_set.append(agent)
        return connected_set


    def find_occupied(self, connected_agents):
        """ 
        Find all the sites that are already occupied on the selected
        binding partner and its connected agents.
        """
    
        occupied = []
        for agent in connected_agents:
            agent_type = self.agents_dictionary[agent]["type"]
            agent_loc = self.agents_dictionary[agent]["loc"]
            sites = self.agents_dictionary[agent]["sites"]
            for site in sites.keys():
                occupied.append([agent, agent_type, agent_loc, site])
        return occupied


    def find_possible(self, connected_agents):
        """ Find all possible sites according to the interaction matrix. """
    
        possible = []
        for agent in connected_agents:
            agent_type = self.agents_dictionary[agent]["type"]
            agent_loc = self.agents_dictionary[agent]["loc"]
            matrix_row = self.protein_ids.index(agent_type)
            inter_vector = self.inter_matrix[matrix_row]
            for i in range(len(inter_vector)):
                inter = inter_vector[i]
                if inter == 1:
                    possible.append([agent, agent_type, agent_loc,
                                     self.protein_ids[i].lower()])
        return possible


    def find_available(self, possible, occupied):
        """ Sites still available for extension. """
    
        available = []
        for site in possible:
            if site not in occupied:
                available.append(site)
        return available


    def select_site(self):
        """ Select one if the available sites for extension. """
    
        self.selected_site = random.choice(self.available_sites)


    def find_mutation_site(self, complx1, complx2):
        """ Find the custom mutation site given for test. """

        parenthesis = self.mutation_site.index("(")
        custom_agent = self.mutation_site[:parenthesis]
        custom_site = self.mutation_site[parenthesis+1:-1]
        custom_type = custom_agent
        custom_loc = 1
        if "-" in custom_agent:
            dash = custom_agent.index("-")
            custom_type = custom_agent[:dash]
            custom_loc = int(custom_agent[dash+1:])
        self.selected_site = [custom_agent, custom_type, custom_loc, custom_site]
        if custom_agent in complx1:
            self.random_partner = 0
            self.selected_partner = self.partners_list[0]
            self.unselected_partner = self.partners_list[1]
        elif custom_agent in complx2:
            self.random_partner = 1
            self.selected_partner = self.partners_list[1]
            self.unselected_partner = self.partners_list[0]
        else:
            raise Exception("Mutation site {} not found"
                             .format(self.mutation_site))


    def next_bond_num(self):
        """ Find next binding number. """
        
        bnd_numbers = []
        for sites in self.agents_dictionary.values():
            for bnd_str in sites["sites"].values():
                if "/" in bnd_str:
                    slash = bnd_str.index("/")
                    bnd_num = int(bnd_str[slash+1:])
                else:
                    if bnd_str != ".":
                        bnd_num = int(bnd_str)
                    else:
                        pass
                bnd_numbers.append(bnd_num)
        self.next_bnd = max(bnd_numbers)+1


    def find_index(self, sel_site):
        """ Find the location of an existing agent in agents_list. """
    
        agent_ind = False
        agent_type = sel_site[1]
        agent_loc = sel_site[2]
        position = 0
        for i in range(len(self.agents_list)):
            agent_str = self.agents_list[i]
            parenthesis = agent_str.index("(")
            agent_id = agent_str[:parenthesis]
            if agent_id == agent_type:
                position += 1
                if position == agent_loc:
                    agent_ind = i
        return agent_ind
    # self.source_index = find_index(self.selected_site)


    def new_agents(self):
        """ Create new lists for the agents of the new rules. """
    
        new_agents_list = []
        for agent in self.agents_list:
            new_agents_list.append(agent)
        return new_agents_list
    # agents_free = new_agents(agents_list)


    def add_new_sites(self, agent_ind, new_agents_free, new_agents_bound):
        """ Add binding site to existing agent. """
    
        old_str = self.agents_list[agent_ind]
        new_free = old_str[:-1]+" {}[.])".format(self.selected_site[3])
        new_agents_free[self.source_index] = new_free
        new_bound = old_str[:-1]+" {}[{}])".format(self.selected_site[3], self.next_bnd)
        new_agents_bound[self.source_index] = new_bound    
        return new_agents_free, new_agents_bound
    # agents_free, agents_bound = add_new_sites(agents_free, agents_bound)


    def add_new_ring(self, ring_ind, rng_site, new_agents_ring):
        """ Add binding sites to form a ring with existing agents. """
    
        old_str_sel = self.agents_list[self.source_index]
        old_str_rng = self.agents_list[ring_ind]
        new_sel = old_str_sel[:-1]+" {}[{}])".format(self.selected_site[3], self.next_bnd)
        new_rng = old_str_rng[:-1]+" {}[{}])".format(rng_site[3], self.next_bnd)
        new_agents_ring[self.source_index] = new_sel
        new_agents_ring[ring_ind] = new_rng
        return new_agents_ring
    # agents_ring = add_new_ring(ring_index, ring_site, agents_ring)


    def add_new_agent(self, new_agents_bound):
        """ 
        Insert new protein into the list of agents.
        If the modified partner is on the left, add the new bound protein
        just before the right binding partner. If the modified partner is
        on the right, add the new protein at the end of the list.
        """

        new_protein = "{}({}[{}])".format(self.selected_site[3].upper(),
                                          self.selected_site[1].lower(),
                                          self.next_bnd)
        if self.random_partner == 0:
            for i in range(len(new_agents_bound)-1, -1, -1):
                agent_str = new_agents_bound[i]
                if "/" in agent_str:
                    right_index = i
                    break
            new_agents_bound.insert(right_index, new_protein)
        elif self.random_partner == 1:
            new_agents_bound.append(new_protein)
        return new_agents_bound
    # agents_bound = add_new_agent(agents_bound)


    def create_new_name(self, name, number):
        """ Create a new name for the new rules. """
    
        alpha_num = number
        if number > 9:
            letters = list(string.ascii_lowercase)
            alpha_num = letters[number-10]            
        new_name = "{}{}".format(name, alpha_num)
        return new_name
    # name_free = create_new_name(rule_name, 1)


    def create_new_rule(self, new_agents, new_name, rate):
        """ Assemble the new rules. """
    
        new_str = ", ".join(new_agents)
        new_rule = "'{}' {} @ {} // {}({})".format(new_name, new_str, rate,
                                                   self.selected_site[0],
                                                   self.selected_site[3])
    
        return new_rule
    # rule_free = create_new_rule(agents_free, name_free)


    def random_select(self):
        """ 
        Randomly choose a mutation site. This first requires to look at
        what complexes are involved on the rule.
        """
    
        self.read_inter_matrix()
        self.find_binding_rules()
        self.select_rule()
        self.rule_name_rate()
        self.agents_connections()
        if self.mutation_site == None:
            self.select_partner()
        else:
            complex1 = self.find_complex(self.partners_list[0])
            complex2 = self.find_complex(self.partners_list[1])
            self.find_mutation_site(complex1, complex2)
        self.agents_complex = self.find_complex(self.selected_partner)
        self.occupied_sites = self.find_occupied(self.agents_complex)
        self.possible_sites = self.find_possible(self.agents_complex)
        self.available_sites = self.find_available(self.possible_sites,
                                                   self.occupied_sites)
        if self.mutation_site == None:
            self.select_site()
        else:
            if self.selected_site not in self.available_sites:
                raise Exception("Provided mutation_site already bound "
                                "in rule from select_line.")
        self.next_bond_num()

 
    def build_new_rules(self):
        """ 
        Create two new refined rules. One where the selected_site is
        free and one where it is bound.
        """

        # Create new agent lists.
        self.source_index = self.find_index(self.selected_site)
        agents_free = self.new_agents()
        agents_bound = self.new_agents()
        agents_free, agents_bound = self.add_new_sites(self.source_index,
            agents_free, agents_bound)
        agents_bound = self.add_new_agent(agents_bound)
        # Create names for the new rules.
        name_free = self.create_new_name(self.rule_name, 1)
        name_bound = self.create_new_name(self.rule_name, 2)
        name_uni_free = self.create_new_name(self.unary_name, 1)
        name_uni_bound = self.create_new_name(self.unary_name, 2)
        # Create the new free and bound rules.
        self.rule_free = self.create_new_rule(agents_free, name_free,
                                              self.rule_rate)
        self.rule_bound = self.create_new_rule(agents_bound, name_bound,
                                               self.rule_rate)
        self.uni_free = self.create_new_rule(agents_free, name_uni_free,
                                             self.unary_rate)
        self.uni_bound = self.create_new_rule(agents_bound, name_uni_bound,
                                              self.unary_rate)


    def build_new_rings(self):
        """
        Deal with possible rings. They can be formed if there already is an
        available partner for the chosen site within the existing complex.
        Include the other side of the bond to be formed.
        """

        self.other_complex = self.find_complex(self.unselected_partner)
        self.occupied_other = self.find_occupied(self.other_complex)
        self.possible_other = self.find_possible(self.other_complex)
        self.available_other = self.find_available(self.possible_other,
                                                   self.occupied_other)
        self.possible_ring = []
        for site in self.available_other:
            self.possible_ring.append(site)
        # If the rule is already a ring closure, self.available_sites and
        # self.available_other will contain the same elements (but in a
        # different order). If the rule is not a ring closure,
        # self.available_sites and self.available_other will be completely
        # disjoint.
        # So, if I find at least one common element between
        # self.available_sites and self.available_other, it means that the
        # rule is already a ring closure.
        self.is_closure = False
        for site in self.available_sites:
            if site not in self.available_other:
                self.possible_ring.append(site)
            else:
                self.is_closure = True
        # No site from the chosen_partner should serve to close a ring.
        # It would form an intramolecular bond.
        auto_sites = []
        for site in self.possible_ring:
            if site[0] == self.selected_site[0]:
                auto_sites.append(site)
        for site in auto_sites:
            self.possible_ring.remove(site)
        # Find all possible rings.
        self.ring_matches = []
        for candidate in self.possible_ring:
            candidate_type = candidate[1]
            if (candidate[1].lower() == self.selected_site[3] and
                    candidate[3] == self.selected_site[1].lower()):
                self.ring_matches.append(candidate)
        # Create a rule for each possible ring.
        rule_num = 3 # rule_free and rule_bound are numbers 1 and 2.
        self.rules_ring = []
        self.uni_rules_ring = []
        for ring_site in self.ring_matches:
            new_closure = False
            if self.is_closure == False and ring_site in self.available_other:
                new_closure = True
            ring_index = self.find_index(ring_site)
            agents_ring = self.new_agents()
            agents_ring = self.add_new_ring(ring_index, ring_site, agents_ring)
            name_ring = self.create_new_name(self.rule_name, rule_num)
            uni_name_ring = self.create_new_name(self.unary_name, rule_num)
            if new_closure == True:
                space = name_ring.rfind(" ")
                name_ring = "{} cls{}".format(name_ring[:space],
                                              name_ring[space:])
                curl_open = self.unary_rate.index("{")
                curl_close = self.unary_rate.index("}")
                rate = self.unary_rate[curl_open+1:curl_close]
            else:
                rate = self.rule_rate
            rule_ring = self.create_new_rule(agents_ring, name_ring,
                                             rate)
            uni_rule_ring = self.create_new_rule(agents_ring, uni_name_ring,
                                             self.unary_rate)
            if new_closure == True:
                uni_rule_ring = "//{}".format(uni_rule_ring)
            rule_num += 1
            self.rules_ring.append(rule_ring)
            self.uni_rules_ring.append(uni_rule_ring)


    def build_new_file(self):
        """ Build the Kappa file containing the new mutated rules. """

        if "/" in self.kappa_file:
            slash = self.kappa_file.rfind("/")
            father = self.kappa_file[slash+1:]
        else:
            father = self.kappa_file
        self.new_file = "// Father: {} \n".format(father)
        self.mutated_lines = []
        for i in range(1, len(self.kappa_obj)):
            line = self.kappa_obj[i]
            if i == self.random_line:
                self.mutated_lines.append(i)
                self.mutated_lines.append(i+1)
                self.new_file += "{}\n".format(self.rule_free)
                self.new_file += "{}\n".format(self.rule_bound)
                mut_line = i+1
                added_lines = 1
                for ring_rule in self.rules_ring:
                    mut_line += 1
                    added_lines += 1
                    self.mutated_lines.append(mut_line)
                    self.new_file += "{}\n".format(ring_rule)
            elif i == self.unary_line:
                self.mutated_lines.append(i+added_lines)
                self.mutated_lines.append(i+added_lines+1)
                if self.unary_commented == True:
                    self.new_file += "//"
                self.new_file += "{}\n".format(self.uni_free)
                if self.unary_commented == True:
                    self.new_file += "//"
                self.new_file += "{}\n".format(self.uni_bound)
                mut_line = i+added_lines+1
                for uni_ring_rule in self.uni_rules_ring:
                    mut_line += 1
                    self.mutated_lines.append(mut_line)
                    if self.unary_commented == True:
                        self.new_file += "//"
                    self.new_file += "{}\n".format(uni_ring_rule)
            else:
                self.new_file += "{}".format(line)


    def __str__(self):
        """
        String representation of the MutateRule class.
        Show the original selected rule and the mutated rules.
        """

        res = ("Selected rule (line {} from file {}) :\n\n"
               .format(self.random_line, self.kappa_file))
        res += "{}\n\n".format(self.selected_rule)
        res += ("New rules with mutations on selected site {}({}) :\n\n"
                .format(self.selected_site[0], self.selected_site[3]))
        res += "{}\n{}\n".format(self.rule_free, self.rule_bound)
        for ring_rule in self.rules_ring:
            res += "{}\n".format(ring_rule)
        return res


    def __repr__(self):
        """ Representation of the MutateRule object. """

        return self.new_file


class MutateRate:
    """
    Change the rate of rules from a Kappa file.
    Either put a random initial value to every rate for initialization or
    change the rate from one random rule for evolutionary algorithm.
    Four different rate values are possible: zero, slow, normal and fast.
    """

    def __init__(self, kappa_file, binary_rates, unary_rates,
                 select_lines=None, allrates=False, unchanged_father=False):
        """ Initialize MutateRate class. """

        #                  zero  slow     normal  fast
        self.rate_values = binary_rates
        self.unary_values = unary_rates
        self.rule_arity = ["binary", "unary"]
        self.kappa_file = kappa_file
        self.select_lines = select_lines # Should be a list.
        self.allrates = allrates
        self.unchanged_father = unchanged_father
        self.kappa_obj = open(self.kappa_file, "r").readlines()
        # Run class methods.
        self.change_rates()
        self.build_new_file()

    
    def find_binding_rules(self):
        """ Read input Kappa file and find the range of the binding rules. """

        self.bind_start = 0
        self.bind_end = 0
        self.unary_start = 0
        self.unary_end = 0
        for i in range(len(self.kappa_obj)):
            line = self.kappa_obj[i]
            if "// Binary binding rules." in line:
                self.bind_start = i + 1
            if "// Unary binding rules." in line:
                self.bind_end = i - 1
                self.unary_start = i + 1
            if "// Unbinding rules." in line:
                self.unary_end = i - 1


    def get_rules_rates(self):
        """ Make a list of the rate """

        self.name_list = []
        self.rule_list = []
        self.rate_list = []
        self.rest_list = []
        for i in range(self.bind_start, self.bind_end+1):
            line = self.kappa_obj[i]
            name_begin = line.index("'")
            name_end = line.rfind("'")
            name = line[name_begin+1:name_end]
            rate_begin = line.rfind("@")
            if "cls" in name:
                rate_end = line.rfind("//")
            else:
                rate_end = line.rfind("{")
            rule = line[:rate_begin+2]
            rate = line[rate_begin+2:rate_end-1]
            rest = line[rate_end:-1]
            self.name_list.append(name)
            self.rule_list.append(rule)
            self.rate_list.append(rate)
            self.rest_list.append(rest)
        #self.unary_name_list = []
        self.unary_rule_list = []
        self.unary_rate_list = []
        self.unary_rest_list = []
        for i in range(self.unary_start, self.unary_end+1):
            line = self.kappa_obj[i]
            rate_begin = line.rfind("@")
            unary_begin = line.rfind("{")
            unary_end = line.rfind("}")
            rule = line[:rate_begin+2]
            rate = line[unary_begin+1:unary_end]
            rest = line[unary_end+1:-1]
            self.unary_rule_list.append(rule)
            self.unary_rate_list.append(rate)
            self.unary_rest_list.append(rest)
        if len(self.rate_list) != len(self.unary_rate_list):
            raise Exception("Number of binary and unary rules "
                            "should be the same.")
        self.all_closures = True
        for name in self.name_list:
            if "cls" not in name:
                self.all_closures = False
                break


    def select_rule(self):
        """ Randomly select a binding rule of chosen_arity. """

        if self.chosen_arity == "binary":
            self.random_line = random.randrange(self.bind_start,
                                                self.bind_end+1)
        elif self.chosen_arity == "unary":
            self.random_line = random.randrange(self.unary_start,
                                                self.unary_end+1)


    def fresh_binary_rates(self):
        """
        Make a fresh copy of the list containing available binary rate values.
        """

        self.available_rates = []
        for val in self.rate_values:
            self.available_rates.append(val)


    def fresh_unary_rates(self):
        """
        Make a fresh copy of the list containing available unary rate values.
        """

        self.available_rates = []
        for val in self.unary_values:
            self.available_rates.append(val)


    def change_rates(self):
        """ Assign new rates values """

        self.find_binding_rules()
        self.get_rules_rates()
        if self.allrates == True:
            self.change_all()
        else:
            if self.select_lines == None:
                self.change_random()
            elif self.select_lines != None:
                self.change_selected()


    def change_all(self):
        """
        Change all rates.
        A new value can be equal to the old one.
        """

        for i in range(len(self.rate_list)):
            self.rule_index = i
            rule_str = self.rule_list[self.rule_index]
            if "cls" in rule_str:
                self.fresh_unary_rates()
            else:
                self.fresh_binary_rates()
            random_rate = random.choice(self.available_rates)
            new_rate = "{}".format(random_rate)
            self.rate_list[i] = new_rate
        for i in range(len(self.unary_rate_list)):
            self.rule_index = i
            rule_str = self.unary_rule_list[self.rule_index]
            if rule_str[:2] != "//":
                corr_bin_rate = float(self.rate_list[self.rule_index])
                if corr_bin_rate == 0.0:
                    random_rate = 0
                else:
                    self.fresh_unary_rates()
                    random_rate = random.choice(self.available_rates)
                new_rate = "{}".format(random_rate)
                self.unary_rate_list[self.rule_index] = new_rate


    def change_random(self):
        """
        Change a single random rate.
        The new value will be different from the old one.
        """

        if self.all_closures == True:     # If all binary rules became ring
            self.chosen_arity = "binary"  # closures, all unary are commented.
        else:
            self.chosen_arity = random.choice(self.rule_arity)
        if self.chosen_arity == "binary":
            self.select_rule()
            self.rule_index = self.random_line - self.bind_start
            rule_str = self.rule_list[self.rule_index]
            if "cls" in rule_str:
                self.fresh_unary_rates()
            else:
                self.fresh_binary_rates()
            self.old_rate = float(self.rate_list[self.rule_index])
            self.available_rates.remove(self.old_rate)
            random_rate = random.choice(self.available_rates)
            new_rate = "{}".format(random_rate)
            self.rate_list[self.rule_index] = new_rate
        elif self.chosen_arity == "unary":
            # Cannot select a commented unary rule.
            rule_found = False
            while rule_found == False:
                self.select_rule()
                self.rule_index = self.random_line - self.unary_start
                rule_str = self.unary_rule_list[self.rule_index]
                if rule_str[:2] != "//":
                    rule_found = True
            self.fresh_unary_rates()
            self.old_rate = float(self.unary_rate_list[self.rule_index])
            self.available_rates.remove(self.old_rate)
            random_rate = random.choice(self.available_rates)
            # We do not allow non zero unary rates if the
            # corresponding binary rate is 0.
            binary_rate = float(self.rate_list[self.rule_index])
            if binary_rate == 0.0 and random_rate != 0.0:
                # Change an other binary rate instead.
                self.chosen_arity = "binary"
                # Cannot select a ring closure at that point.
                rule_found = False
                while rule_found == False:
                    self.select_rule()
                    self.rule_index = self.random_line - self.bind_start
                    name_str = self.name_list[self.rule_index]
                    if "cls" not in name_str:
                        rule_found = True
                self.fresh_binary_rates()
                self.rule_index = self.random_line - self.bind_start
                self.old_rate = float(self.rate_list[self.rule_index])
                self.available_rates.remove(self.old_rate)
                random_rate = random.choice(self.available_rates)
                new_rate = "{}".format(random_rate)
                self.rate_list[self.rule_index] = new_rate
            else:
                new_rate = "{}".format(random_rate)
                self.unary_rate_list[self.rule_index] = new_rate


    def change_selected(self):
        """
        Change the rate of each given line.
        The new value can be equal to the old value.
        Any unimolecular rate is set to 0 if the corresponding
        bimolecular rate is 0.
        """

        self.modified_lines = []
        for line in self.select_lines:
            if line >= self.bind_start and line <= self.bind_end:
                self.rule_index = line - self.bind_start
                #self.old_rate = float(self.rate_list[self.rule_index])
                if "cls" in self.name_list[self.rule_index]:
                    self.fresh_unary_rates()
                else:
                    self.fresh_binary_rates()
                random_rate = random.choice(self.available_rates)
                new_rate = "{}".format(random_rate)
                self.rate_list[self.rule_index] = new_rate
                self.modified_lines.append(line)
            if line >= self.unary_start and line <= self.unary_end:
                self.rule_index = line - self.unary_start
                #self.old_rate = float(self.rate_list[self.rule_index])
                rule_str = self.unary_rule_list[self.rule_index]
                if rule_str[:2] != "//":
                    corr_bin_rate = float(self.rate_list[self.rule_index])
                    if corr_bin_rate == 0.0:
                        random_rate = 0
                    else:
                        self.fresh_unary_rates()
                        random_rate = random.choice(self.available_rates)
                    new_rate = "{}".format(random_rate)
                    self.unary_rate_list[self.rule_index] = new_rate
                    self.modified_lines.append(line)
        self.left_out = []
        for line in self.select_lines:
            if line not in self.modified_lines:
                self.left_out.append(line)


    def build_new_file(self):
        """ Build the Kappa file containing the new changed rates. """

        if self.unchanged_father == False:
            if "/" in self.kappa_file:
                slash = self.kappa_file.rfind("/")
                father = self.kappa_file[slash+1:]
            else:
                father = self.kappa_file
            self.new_file = "// Father: {} \n".format(father)
            first_line = 1
        elif self.unchanged_father == True:
            self.new_file = ""
            first_line = 0
        for i in range(first_line, len(self.kappa_obj)):
            line = self.kappa_obj[i]
            if i >= self.bind_start and i <= self.bind_end:
                j = i - self.bind_start
                self.new_file += ("{}{} {}\n"
                    .format(self.rule_list[j], self.rate_list[j],
                            self.rest_list[j]))
            elif i >= self.unary_start and i <= self.unary_end:
                j = i - self.unary_start
                self.new_file += ("{}0 {{{}}}{}\n"
                    .format(self.unary_rule_list[j], self.unary_rate_list[j],
                            self.unary_rest_list[j]))
            else:
                self.new_file += "{}".format(line)


    def __str__(self):
        """ String representation of the MutateRate class. """
        
        res = ""
        if self.allrates == True:
            res += ("Randomly modified all rates from file {}\n"
                    .format(self.kappa_file))
        elif self.allrates == False:
            if self.select_lines == None:
                res += ("Modified rate of rule from line {} of file {}\n"
                        .format(self.random_line, self.kappa_file))
                if self.chosen_arity == "binary":
                    res += ("{}{} --> {}\n"
                            .format(self.rule_list[self.rule_index],
                                    self.old_rate,
                                    self.rate_list[self.rule_index]))
                elif self.chosen_arity == "unary":
                    res += ("{}{} --> {}\n"
                            .format(self.unary_rule_list[self.rule_index],
                                    self.old_rate,
                                    self.unary_rate_list[self.rule_index]))
            elif self.select_lines != None:
                res += ("Modified rate of rule from lines {} of file {}\n"
                        .format(self.modified_lines, self.kappa_file))
                res += ("Lines {} were left untouched since they are commented"
                        .format(self.left_out))
        return res


    def __repr__(self):
        """ Representation of the MutateRate object. """

        return self.new_file

