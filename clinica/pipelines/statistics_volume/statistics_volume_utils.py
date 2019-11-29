# coding: utf8
import pandas as pds
from clinica.utils.exceptions import ClinicaException


def get_group_1_and_2(file_list, csv_file, contrast):

    csv = pds.read_csv(csv_file, sep='\t')
    columns = list(csv.columns)

    if contrast not in columns:
        raise ClinicaException(contrast + ' is not present in ' + csv_file)
