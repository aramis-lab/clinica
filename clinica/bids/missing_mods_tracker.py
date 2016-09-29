
class MissingModsTracker:
    """
    Class used for tracking the number of missing modalities in a database
    """
    mods_missing = {
            'DTI': 0,
            'fMRI': 0,
            'Map': 0,
            'MapPh': 0,
            'FLAIR': 0,
            'T1': 0
        }
    mods_used = []

    def add_missing_mod(self, mod):
        """
        Increase the number of missing files for the input modality.

        Args:
            mod: modality missing
        """
        self.mods_missing[mod] = self.mods_missing[mod]+1

        if mod not in self.mods_used:
            self.mods_used.append(mod)

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
        return self.mods_missing







