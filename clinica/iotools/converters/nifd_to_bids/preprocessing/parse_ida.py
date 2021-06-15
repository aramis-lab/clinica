# coding: utf8


import contextlib
import os


def csv_to_tsv(file, path_out):
    import os

    fich = open(file, "r")
    out = open(os.path.join(path_out, "ida.tsv"), "w")
    s = ""
    for line in fich.readlines():
        line_mod = ""
        line = line.split(",")
        for elt in line:
            line_mod += elt + "\t"
        line_mod = line_mod[:-1]
        s += line_mod

    out.write(s)
    out.close()
    fich.close()


def first_process(path_out):
    import os

    path_ida = os.path.join(path_out, "ida.tsv")

    def get_transf_date(pat_name):
        dates = []
        for elt in pat_dic[pat]:
            line = elt.split("\t")
            date = line[5]
            if date not in dates:
                dates.append(date)
            val = float("inf")
            for date in dates:
                month = int(date.split("/")[0]) + 12 * int(date.split("/")[2])
                if month < val:
                    low_date = date
                    val = month

            transf_date = dict()
            for date in dates:
                transf_date[date] = (
                    int(date.split("/")[0]) + 12 * int(date.split("/")[2]) - val
                )
        return transf_date

    def modif_line(pat_name, final_sol):
        transf_date = get_transf_date(pat_name)
        for line in pat_dic[pat]:
            line = line.split("\t")
            line[4] = "Month " + str(transf_date[line[5]])

            for elt in line:
                final_sol += elt + "\t"
            final_sol = final_sol[:-1]
        return final_sol

    fich = open(path_ida, "r")
    final_sol = fich.readline()
    patients = []
    pat_dic = dict()
    for line in fich.readlines():
        pat = line.split("\t")[0]
        if pat not in patients:
            patients.append(pat)
            pat_dic[pat] = [str(line)]
        else:
            pat_dic[pat].extend([str(line)])

    for pat in patients:
        final_sol = modif_line(pat, final_sol)

    fich.close()
    fich = open(path_ida, "w")
    fich.write(final_sol)
    fich.close()


def second_process(path_out):
    import os

    import pandas as pd

    path_ida = os.path.join(path_out, "ida.tsv")

    dfIda = pd.read_csv(path_ida, sep="\t")
    pat_list = list(set(list(dfIda["Subject ID"])))
    pat_list.sort()

    dfFinal = dfIda.iloc[0:0]

    for pat in pat_list:

        dfPat = dfIda[dfIda["Subject ID"] == pat]
        out1 = list(set(list(dfPat["Visit"])))
        out1 = [int(month.split(" ")[-1]) for month in out1]
        out1.sort()

        if pat not in ["1_S_0067", "1_S_0119", "1_S_0091", "1_S_0093", "2_S_0010"]:
            map = [i - i % 6 if (i % 6 < 6 - i % 6) else i + 6 - i % 6 for i in out1]
        else:
            map = [i - i % 6 if (i % 6 <= 6 - i % 6) else i + 6 - i % 6 for i in out1]

        if pat == "1_S_0030":
            map = [0, 6, 12]
        if pat == "1_S_0316":
            map = [0, 6, 12]
        if pat == "1_S_0334":
            map = [0, 12, 18, 24]
        if pat == "1_S_0349":
            map = [0, 6]
        if pat == "1_S_0354":
            map = [0, 6]

        # Prints the cases where the mapping process creates an overlap (month 0, month 1 will be both be mapped to month 0 for instance)
        # if len(map) != len(set(map)):
        #     print(out1)
        #     print(map)
        #     print('-------------')

        # Prints the cases where an overlap creates a conflict
        # if pat in ['1_S_0030', '1_S_0316', '1_S_0334', '1_S_0349', '1_S_0354']:
        #     print(out1)
        #     print(map)
        #     print('----------')

        dfPat["Visit"] = dfPat["Visit"].apply(
            lambda x: "Month " + str(map[out1.index(int(x.split(" ")[-1]))])
        )
        dfFinal = dfFinal.append(dfPat)

    def write(df, path, name):
        import os

        name = os.path.join(path, name) + ".tsv"
        df = df.fillna("n/a")
        df.to_csv(sep="\t", path_or_buf=name, index=False)

    write(dfFinal, path_out, "ida")


def supress_stdout(func):
    def wrapper(*a, **ka):
        with open(os.devnull, "w") as devnull:
            with contextlib.redirect_stdout(devnull):
                func(*a, **ka)

    return wrapper


@supress_stdout
def process_ida(path_ida, path_out=None):
    """
    Builds 'ida.tsv' from 'idaSearch_all.csv'

    Args:
        path_ida: Path to the original ida file 'idaSearch_all.csv'
        path_out: Path to the output directory for 'ida.tsv', path_out is None, the ouput directory will be the same as the orginal ida file
    """
    import os

    if path_out is None:
        path_out = path_ida.split("/")[:-1]
        path_out = "/" + os.path.join(*path_out)

    # Process from original csv to tsv
    csv_to_tsv(path_ida, path_out)

    # First processing on ida.tsv, recomputes the 'Visit' column to avoid overlaps in IRM visits
    first_process(path_out)

    # Second processing on ida.tsv, maps the 'Visit' column to multiples of 6, still avoiding overlaps
    second_process(path_out)
