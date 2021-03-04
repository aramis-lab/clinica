# coding: utf8


class Descriptor(object):
    """
    Class that handles the conversion of a raw image name from the NIFD dataset,
        handles the BIDS format and the hierarchy of quality
        (T1 images with DIS3D will be selected over basic T1 images,
        the descriptor in charge of finding T1 images with DIS3D will have a higher priority value than the descriptor in charge of finding images with just T1 in its name)
    A descriptor object is created from a dictionnary written in the config_dcm2bids.json file
    To add a descriptor to the descriptors list that is used in convert_images, add it to the config_dcm2bids.json file
    """

    def __init__(self, dic_des):
        self.dataType = None
        self.modalityLabel = None
        self.customLabels = None
        self.priority = 0

        self.criteria = False
        self.Modality = None
        self.SeriesDescription = None

        if "dataType" in dic_des:
            self.dataType = dic_des["dataType"]
        if "modalityLabel" in dic_des:
            self.modalityLabel = dic_des["modalityLabel"]
        if "customLabels" in dic_des:
            self.customLabels = dic_des["customLabels"]
        if "priority" in dic_des:
            self.priority = int(dic_des["priority"])

        # Criteria attributes
        if "criteria" in dic_des:
            self.criteria = True
            if "Modality" in dic_des["criteria"]:
                self.Modality = dic_des["criteria"]["Modality"]
            if "SeriesDescription" in dic_des["criteria"]:
                self.SeriesDescription = dic_des["criteria"]["SeriesDescription"]

    def describes(self, str_image):
        """
        Uses self.SeriesDescription to determine whether self describes str_image
        For instance, if self.SeriesDescription = "*FLAIR*^*DIS3D:*flair*^*DIS3D"
        every image containing the string "FLAIR" and ending with "DIS3D" will make this method return True
        an image containing the string "flair" and ending with "DIS3D" will also return True
        '*' stands for any string
        '^' stands for AND
        ':' stands for OR

        Args:
            str_image: Image name

        Returns:
            Boolean: True if str_image corresponds to the pattern written in self.SeriesDescription, else False
        """

        def test_desc(str_image, desc):
            if desc[0] == "*":
                if desc[-1] == "*":
                    return desc[1:-1] in str_image
                else:
                    return str_image[-len(desc) + 1 :] == desc[1:]

            else:
                if desc[-1] == "*":
                    return str_image[: len(desc) - 1] == desc[:-1]
                else:
                    return str_image == desc

        descriptions = self.SeriesDescription.split(":")
        for desc in descriptions:
            sub_desc = desc.split("^")
            bool = True
            for sub in sub_desc:
                bool = bool and test_desc(str_image, sub)

            if bool:
                return True
        return False

    def get_bids_info(self):
        s = ""
        if self.customLabels is not None:
            s += self.customLabels + "_"
        if self.modalityLabel is not None:
            s += self.modalityLabel + "_"

        if len(s) == 0:
            return s
        return s[:-1]

    def __str__(self):
        s = ""

        if self.dataType is not None:
            s += "dataType : " + self.dataType + ", "
        if self.modalityLabel is not None:
            s += "modalityLabel : " + self.modalityLabel + ", "
        if self.customLabels is not None:
            s += "customLabels : " + self.customLabels + ", "
        if self.priority != 0:
            s += "priority : " + str(self.priority) + ", "

        # Criteria attributes
        if self.criteria:
            s += "\n(criteria) "
            if self.Modality is not None:
                s += "Modality : " + self.Modality + ", "
            if self.SeriesDescription is not None:
                s += "SeriesDescription : " + self.SeriesDescription + ", "

        if len(s) > 2:
            return s[:-2]
        return s
