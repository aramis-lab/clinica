import pandas

def print_statistics(summary_file, num_subjs, ses_aval, mmt):
    missing_list = mmt.get_missing_list()
    summary_file.write('Number of subjects converted: ' + str(num_subjs) + '\n')
    summary_file.write('Sessions available: '+ses_aval[0] +' '+ses_aval[1]+'\n')


    for ses in ses_aval:
        ses_miss = missing_list[ses]['session']
        ses_found = num_subjs - ses_miss
        perc_ses_found = ses_found*100/num_subjs
        summary_file.write('Number of sessions '+ses+" found: "+str(ses_found)+' ('+str(perc_ses_found)+'%)\n')

    summary_file.write('****************************************\n')
    summary_file.write('Number of missing modalities for each session (more details inside converter.log):\n')


    for ses in ses_aval:
         summary_file.write('\n-'+ses+'-\n')
         for mod in missing_list[ses]:
             if mod != 'session':
                 num_miss_mod = missing_list[ses][mod][0]
                 percentage_missing = round((num_miss_mod*100/float(num_subjs - missing_list[ses]['session'])),2)
                 summary_file.write(mod+': '+str(num_miss_mod)+' ('+str(percentage_missing)+'%) \n')


    summary_file.write('\n\nNumber of incomplete modalities for each session (more details inside converter.log):\n')
    for ses in ses_aval:
         summary_file.write('\n-'+ses+'-\n')
         for mod in missing_list[ses]:
             if mod != 'session':
                 num_miss_mod = missing_list[ses][mod][1]
                 if num_miss_mod !=0:
                    percentage_missing = round((num_miss_mod*100/float(num_subjs - missing_list[ses]['session'])),2)
                 else:
                    percentage_missing = 0
                 summary_file.write(mod+': '+str(num_miss_mod)+' ('+str(percentage_missing)+'%) \n')


class MissingModsTracker:
    """
    Class used for tracking the number of missing modalities in a database
    """
    def __init__(self, ses):
        self.missing = {}
        self.mods_missing = {
            'DTI': 0,
            'fMRI': 0,
            'Fieldmap': 0,
            'FLAIR': 0,
            'T1': 0
        }
        # Different session
        if ses != '':
            self.ses= ses
            for s in ses:
                self.missing.update({s : {'session': 0,
                                          'DTI': [0,0],
                                          'fMRI': [0,0],
                                          'Fieldmap': [0,0],
                                          'FLAIR': [0,0],
                                          'T1': [0,0]}
                                })
        self.mods_used = []

    def add_missing_mod(self, mod, ses):
        """
        Increase the number of missing files for the input modality.

        Args:
            mod: modality missing
        """
        self.mods_missing[mod] = self.mods_missing[mod]+1
        (self.missing[ses][mod])[0] = (self.missing[ses][mod][0]) +1

        if mod not in self.mods_used:
            self.mods_used.append(mod)

    def increase_missing_ses(self, ses):
        self.missing[ses] = self.missing[ses]+1


    def add_incomplete_mod(self, mod, ses):
        (self.missing[ses][mod])[1] = (self.missing[ses][mod][1]) + 1


    def incr_missing_session(self, ses):
        self.missing[ses]['session'] +=1

    def remove_mods_unused(self):
        """
        Remove all the modalities not included into the list mods_used
        """
        keys = self.mods_missing.keys()
        for key in keys:
            if key not in self.mods_used:
                del self.mods_missing[key]

    def get_missing_list(self):
        """
        Return the hash map mods_missing.

        Returns:
             The hash map containing the list of missing files.

        """
        return self.missing









