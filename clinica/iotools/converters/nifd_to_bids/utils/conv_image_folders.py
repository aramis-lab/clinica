# coding: utf8


from clinica.iotools.converters.nifd_to_bids.utils.descriptor import Descriptor


class Old_image_folder(object):
    def __init__(self, strg):

        self.name = strg
        self.new_name = False
        self.old_name = None

        self.bids_info = None

    def get_new_name(self, desc_list):
        """Find the descriptor that describes the current image with the highest priority.

        Args:
            desc_list: list of descriptor instances

        Returns:
            res_desc: Single descriptor instance
        """
        pos_desc = [desc for desc in desc_list if desc.describes(self.name)]
        if len(pos_desc) == 0:
            return None
        else:
            max_priority = 0
            res_desc = None
            for desc in pos_desc:
                if desc.priority > max_priority:
                    max_priority = desc.priority
                    res_desc = desc
            self.old_name = self.name
            self.name = str(res_desc)
            self.new_name = True
            self.bids_info = res_desc.get_bids_info()
            return res_desc


def get_all_med_name(path_dataset):
    """Create a list 'medical_images' containing all possible medical images' names in the NIFD dataset.

    Args:
        path_dataset: path to the NIFD dataset

    Returns:
        medical_images = ['t1_mprage', 't1_mprage_S3_DIS3D', ...]
    """
    import os

    subfolders = [f.path for f in os.scandir(path_dataset) if f.is_dir()]
    medical_images = []
    for sub in subfolders:
        subfolders2 = [f.path for f in os.scandir(sub) if f.is_dir()]
        for sub2 in subfolders2:
            medical_image_name = os.path.basename(sub2)
            if medical_image_name not in medical_images:
                medical_images.append(medical_image_name)
    return medical_images


def get_descriptors(path_root):
    """Load the descriptors.

    Args:
        path_root: path to the folder containing the 'config_dcm2bids.json' file

    Returns:
        descriptors: list of descriptor instances
    """
    import json
    import os

    with open(os.path.join(path_root, "config_dcm2bids.json")) as f:
        data = json.load(f)

    descriptors = []
    for i in data["descriptions"]:
        des = Descriptor(i)
        descriptors.append(des)

    return descriptors


def dict_conversion(medical_images, descriptors):
    """
    For each medical image in the dataset, find whether its name is described by one of the descriptor object.
    If it is, the descriptor with the highest priority is selected.
    The medical image is then added to the dictionnary as a key with its associated descriptor and its modality label.

    Args:
        medical_images: list containing every medical image name
        descriptors: list of descriptor instances

    Returns:
        equivalences: dictionary, equivalences['medical_image_name'] = (Descriptor_instance, modalityLabel)
    """
    equivalences = dict()

    for med_name in medical_images:
        old_name = Old_image_folder(med_name)
        desc = old_name.get_new_name(descriptors)
        if old_name.new_name:
            equivalences[old_name.old_name] = (desc, old_name.bids_info)

    return equivalences
