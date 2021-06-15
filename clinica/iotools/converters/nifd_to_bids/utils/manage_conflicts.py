# coding: utf8


class Manage_conflicts:
    def __init__(self, path_choices):

        self.dic = self.parse_choices(path_choices)

    def parse_choices(self, path_choices):
        f = open(path_choices, "r")
        dic = dict()
        for line in f.readlines():
            line = line.split("-")
            dic[line[0]] = line[1]

        return dic

    def make_decision(self, list_med_names):
        list_med_names.sort()
        return self.dic[str(list_med_names)]

    def __str__(self):
        s = ""
        for key in self.dic:
            s += str(key) + " -> " + str(self.dic[key])

        return s
