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
    def get_name_atlas(self):pass
    """
    Returns the name of the atlas (as defined in BIDS/CAPS specifications).
    """

    def get_spatial_resolution(self):
        import nibabel as nib

        img_map = nib.load(self.get_atlas_map())
        img_labels = nib.load(self.get_atlas_labels())
        voxels_map = img_map.header.get_zooms()
        voxels_labels = img_labels.header.get_zooms()
        if voxels_map != voxels_labels:
            raise Exception(
                "Spatial resolution of labels and map image from %s atlas mismatch" % (self.get_name_atlas()))
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
    """
    Returns the spatial resolution of the atlas (in format "XxXxX" e.g. 1x1x1 or 1.5x1.5x1.5).
    """

    @abc.abstractmethod
    def get_atlas_labels(self):pass
    """
    Returns the image with the different ROIs.
    """

    @abc.abstractmethod
    def get_atlas_map(self):pass
    """
    Returns the map associated to the atlas (e.g. T1, FA map from DTI, etc.).
    """

    @abc.abstractmethod
    def get_roi_name(self):pass
    """
    Returns the spatial resolution of the atlas (in format "XxXxX" e.g. 1x1x1 or 1.5x1.5x1.5).
    """




class Atlas_JHUDTI81_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUDTI81"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-labels-2mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUDTI81_ROI.tsv')




class Atlas_JHUDTI81_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUDTI81"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-labels-1mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUDTI81_ROI.tsv')




class Atlas_JHUTracts0_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUTracts0"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class Atlas_JHUTracts0_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUTracts0"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr0-2mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')


class Atlas_JHUTracts25_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUTracts25"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class Atlas_JHUTracts25_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUTracts25"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr25-2mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class Atlas_JHUTracts50_1mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUTracts50"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr50-1mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-1mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')



class Atlas_JHUTracts50_2mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    def get_name_atlas(self): return "JHUTracts50"

    def get_atlas_labels(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-tracts-maxprob-thr50-2mm.nii.gz')

    def get_atlas_map(self):
        import os
        FSLDIR = os.environ.get('FSLDIR', '')
        if not FSLDIR:
            raise Exception('FSLDIR variable from FSL software is not set')
        return os.path.join(FSLDIR, 'data', 'atlases', 'JHU', 'JHU-ICBM-FA-2mm.nii.gz')

    def get_roi_name(self):
        import os
        CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
        if not CLINICA_HOME:
            raise Exception('CLINICA_HOME variable from Clinica software is not set')
        return os.path.join(CLINICA_HOME, 'clinica', 'resources', 'atlases', 'JHUTract_ROI.tsv')





class AtlasLoader:
    def __init__(self, atlases=None):
        self.atlas = {}
        if atlases:
            for atlas in atlases: self.add_atlas(atlas)

    def add_atlas(self,atlas):
        if not isinstance(atlas, AtlasAbstract):
            raise Exception("Atlas element must be an AtlasAbstract type")

    def get_atlases(self): return self.atlas
    
    

