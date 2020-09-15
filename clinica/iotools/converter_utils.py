# coding: utf8


def print_statistics(summary_file, num_subjs, ses_aval, mmt):
    """
    Print to a given input file statistics about missing files and modalities
    in a dataset.  This metod takes in input a MissingModsTracker object (mmt)
    that contains the number of missing modalities for each session and the
    number of missing sessions for each subject.

    Args:
        summary_file: path of the output file where write.
        num_subjs: number of subjects.
        ses_aval: list of sessions available.
        mmt: object MissingModsTracker
    """
    missing_list = mmt.get_missing_list()
    summary_file.write('**********************************************\n')
    summary_file.write('Number of subjects converted: ' + str(num_subjs) + '\n')
    summary_file.write('Sessions available: '+ses_aval[0] + '\n')

    for ses in ses_aval:
        ses_miss = missing_list[ses]['session']
        ses_found = num_subjs - ses_miss
        perc_ses_found = ses_found*100/num_subjs
        summary_file.write('Number of sessions '+ses+" found: "+str(ses_found)+' ('+str(perc_ses_found)+'%)\n')

    summary_file.write('**********************************************\n\n')
    summary_file.write('Number of missing modalities for each session:\n')

    for ses in ses_aval:
        summary_file.write('\n' + ses + '\n')
        for mod in missing_list[ses]:
            if mod != 'session':
                num_miss_mod = missing_list[ses][mod]
                percentage_missing = round((num_miss_mod*100/float(num_subjs - missing_list[ses]['session'])), 2)
                summary_file.write(mod+': ' + str(num_miss_mod) + ' ('+str(percentage_missing) + '%) \n')


def has_one_index(index_list):
    if len(index_list) == 1:
        return index_list[0]
    if len(index_list) == 0:
        return -1
    if len(index_list) > 1:
        raise('Multiple indexes found')


class MissingModsTracker:
    """
    Class used for tracking the number of missing modalities in a database
    """
    def __init__(self, ses, mod_list=False):
        self.missing = {}
        self.ses = ses
        if mod_list:
            for s in ses:
                self.missing.update({s: {'session': 0}})
                for mod in mod_list:
                    self.missing[s].update({mod: 0})
        else:
            for s in ses:
                self.missing.update({
                    s: {'session': 0,
                        'dwi': 0,
                        'func': 0,
                        'fieldmap': 0,
                        'flair': 0,
                        't1w': 0}
                    })

    def add_missing_mod(self, ses, mod):
        """
        Increase the number of missing files for the input modality.

        Args:
            ses: name of the session
            mod: modality missing
        """
        self.missing[ses][mod] += 1

    def increase_missing_ses(self, ses):
        self.missing[ses]['session'] += 1

    def get_missing_list(self):
        """
        Return the hash map mods_missing.

        Returns:
             The hash map containing the list of missing files.

        """
        return self.missing
