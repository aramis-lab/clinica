#!/usr/bin/env python
# -*- coding: utf-8 -*-
import abc


class AtlasAbstract:
    """
    Abstract class for Atlas handling.

    Naming convention for children classes of AtlasAbstract: <name_atlas>[_<resolution>][_<map>]
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_name_atlas(self):
        """
        Returns the name of the atlas (as defined in BIDS/CAPS specifications).
        """
        pass

    def get_spatial_resolution(self):
        """
        Returns the spatial resolution of the atlas (in format "XxXxX" e.g. 1x1x1 or 1.5x1.5x1.5).
        """
        import nibabel as nib

        img_map = nib.load(self.get_atlas_map())
        img_labels = nib.load(self.get_atlas_labels())
        voxels_map = img_map.header.get_zooms()
        voxels_labels = img_labels.header.get_zooms()
        if (voxels_map[0] != voxels_labels[0]) or \
                (voxels_map[1] != voxels_labels[1]) or \
                (voxels_map[2] != voxels_labels[2]):
        #if voxels_map != voxels_labels:
            print "Spatial resolution of labels and map image from %s atlas mismatch" % (self.get_name_atlas())
            #raise Exception(
             #       "Spatial resolution of labels and map image from %s atlas mismatch" % (self.get_name_atlas()))
        #else:
        # Will display integers without decimals
        if int(voxels_map[0]) == voxels_map[0]:
            s_x = str(int(voxels_map[0]))
        else:
            s_x = str(voxels_map[0])
        if int(voxels_map[1]) == voxels_map[1]:
            s_y = str(int(voxels_map[1]))
        else:
            s_y = str(voxels_map[1])
        if int(voxels_map[2]) == voxels_map[2]:
            s_z = str(int(voxels_map[2]))
        else:
            s_z = str(voxels_map[2])

        return s_x + "x" + s_y + "x" + s_z

    @abc.abstractmethod
    def get_atlas_labels(self):
        """
        Returns the image with the different labels/ROIs.
        """
        pass

    @abc.abstractmethod
    def get_atlas_map(self):
        """
        Returns the map associated to the atlas (e.g. T1, FA map from DTI, etc.).
        """
        pass

    @abc.abstractmethod
    def get_tsv_roi(self):
        """
        Returns the TSV file containing the ROI (regions of interest) of the atlas.
        """
        pass

    def get_index(self):
        import nibabel as nib
        import numpy as np
        img_labels = nib.load(self.get_atlas_labels())
        img_labels = img_labels.get_data()
        labels = list(set(img_labels.ravel()))
        index_vector=np.zeros(len(labels))
        for index, n in enumerate(labels):
            index_vector[index] = index
        return index_vector


class JHUDTI81_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUDTI81"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-labels-2mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUDTI81_ROI.tsv')


class JHUDTI81_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUDTI81"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-labels-1mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUDTI81_ROI.tsv')


class JHUTracts0_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUTracts0"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class JHUTracts0_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUTracts0"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr0-2mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')


class JHUTracts25_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUTracts25"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class JHUTracts25_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUTracts25"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr25-2mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class JHUTracts50_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUTracts50"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr50-1mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')


class JHUTracts50_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "JHUTracts50"

    @staticmethod
    def get_atlas_labels():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr50-2mm.nii.gz')

    @staticmethod
    def get_atlas_map():
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')
    
    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')


class AAL2(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "AAL2"

    @staticmethod
    def get_atlas_labels():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','AAL2.nii')

    @staticmethod
    def get_atlas_map():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'Template_MNI152.nii')

    @staticmethod
    def get_tsv_roi():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'AAL2_ROI.tsv')


class Hammers(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "Hammers"

    @staticmethod
    def get_atlas_labels():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica',  'resources', 'atlases_spm','Hammers.nii')

    @staticmethod
    def get_atlas_map():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','Template_MNI152.nii')

    @staticmethod
    def get_tsv_roi():
        import os.path
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'Hammers_ROI.tsv')


class LPBA40(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "LPBA40"

    @staticmethod
    def get_atlas_labels():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica','resources', 'atlases_spm','LPBA40.nii')

    @staticmethod
    def get_atlas_map():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','Template_MNI152.nii')

    @staticmethod
    def get_tsv_roi():
        import os.path
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'LPBA40_ROI.tsv')


class AICHA(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "AICHA"

    @staticmethod
    def get_atlas_labels():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','AICHA.nii')

    @staticmethod
    def get_atlas_map():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','Template_MNI152.nii')

    @staticmethod
    def get_tsv_roi():
        import os.path
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'AICHA_ROI.tsv')


class Neuromorphometrics(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "Neuromorphometrics"

    @staticmethod
    def get_atlas_labels():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica','resources', 'atlases_spm','Neuromorphometrics.nii')

    @staticmethod
    def get_atlas_map():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','Template_MNI152.nii')

    @staticmethod
    def get_tsv_roi():
        import os.path
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'Neuromorphometrics_ROI.tsv')


class MCALT_ADIR122(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas(): return "MCALT_ADIR122"

    @staticmethod
    def get_atlas_labels():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica','resources', 'atlases_spm','MCALT_ADIR122.nii')

    @staticmethod
    def get_atlas_map():
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm','Template_MNI152.nii')

    @staticmethod
    def get_tsv_roi():
        import os.path
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases_spm', 'MCALT_ADIR122_ROI.tsv')


class AtlasLoader:
    def __init__(self, atlases=None):
        self.atlas = {}
        if atlases:
            for atlas in atlases:
                self.add_atlas(atlas)

    def add_atlas(self, atlas):
        if not isinstance(atlas, AtlasAbstract):
            raise Exception("Atlas element must be an AtlasAbstract type")

    def get_atlases(self): return self.atlas
