__author__ = "Adam Wild"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Adam Wild"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Adam Wild"
__email__ = "adam.wild@icm-institute.org"
__status__ = "Development"

class Manage_conflicts():

    def __init__(self, path_choices):

        self.dic = self.parse_choices(path_choices)

    def parse_choices(self, path_choices):
        f = open(path_choices, 'r')
        dic = dict()
        for line in f.readlines():
            line = line.split('-')
            dic[line[0]] = line[1]

        return dic

    def make_decision(self, list_med_names):
        list_med_names.sort()
        return self.dic[str(list_med_names)]

    def __str__(self):
        s = ''
        for key in self.dic:
            s += str(key) + ' -> ' + str(self.dic[key])

        return s

if __name__=='__main__':
    f = Manage_conflicts('/Users/adam.wild/Desktop/neo_parse_NIFD/unique_conflicts_demo')
    sol = f.make_decision(['T1_mprage_S3_DIS3D', 'T1_mprage_short_S7_DIS3D'])
    print(sol)
