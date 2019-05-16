# coding: utf8


import abc
import os.path as path

import numpy as np
from pandas.io import parsers

from clinica.pipelines.machine_learning import base
import clinica.pipelines.machine_learning.voxel_based_io as vbio
import clinica.pipelines.machine_learning.vertex_based_io as vtxbio
import clinica.pipelines.machine_learning.region_based_io as rbio
import clinica.pipelines.machine_learning.tsv_based_io as tbio
import clinica.pipelines.machine_learning.ml_utils as utils


__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez", "Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class CAPSInput(base.MLInput):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, precomputed_kernel=None):
        """

        Args:
            caps_directory:
            subjects_visits_tsv:
            diagnoses_tsv:
            group_id:
            image_type: 'T1', 'fdg', 'av45', 'pib' or 'flute'
            precomputed_kernel:
        """

        self._caps_directory = caps_directory
        self._group_id = group_id
        self._image_type = image_type
        self._images = None
        self._x = None
        self._y = None
        self._kernel = None

        subjects_visits = parsers.read_csv(subjects_visits_tsv, sep='\t')
        if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
            raise Exception('Subjects and visits file is not in the correct format.')
        self._subjects = list(subjects_visits.participant_id)
        self._sessions = list(subjects_visits.session_id)

        diagnoses = parsers.read_csv(diagnoses_tsv, sep='\t')
        if 'diagnosis' not in list(diagnoses.columns.values):
            raise Exception('Diagnoses file is not in the correct format.')
        self._diagnoses = list(diagnoses.diagnosis)

        if image_type not in ['T1', 'fdg', 'av45', 'pib', 'flute', 'dwi']:
            raise Exception("Incorrect image type. It must be one of the values 'T1', 'fdg', 'av45', 'pib', 'flute' or 'dwi'")

        if precomputed_kernel is not None:
            if type(precomputed_kernel) == np.ndarray:
                if precomputed_kernel.shape == (len(self._subjects), len(self._subjects)):
                    self._kernel = precomputed_kernel
                else:
                    raise Exception("""Precomputed kernel provided is not in the correct format.
                    It must be a numpy.ndarray object with number of rows and columns equal to the number of subjects,
                    or a filename to a numpy txt file containing an object with the described format.""")
            elif type(precomputed_kernel == str):
                self._kernel = np.loadtxt(precomputed_kernel)
            else:
                raise Exception("""Precomputed kernel provided is not in the correct format.
                It must be a numpy.ndarray object with number of rows and columns equal to the number of subjects,
                or a filename to a numpy txt file containing an object with the described format.""")

    @abc.abstractmethod
    def get_images(self):
        """

        Returns: a list of filenames

        """
        pass

    @abc.abstractmethod
    def get_x(self):
        """

        Returns: a numpy 2d-array.

        """
        pass

    def get_y(self):
        """

        Returns: a list of integers. Each integer represents a class.

        """
        if self._y is not None:
            return self._y

        unique = list(set(self._diagnoses))
        self._y = np.array([unique.index(x) for x in self._diagnoses])
        return self._y

    def get_kernel(self, kernel_function=utils.gram_matrix_linear, recompute_if_exists=False):
        """

        Returns: a numpy 2d-array.

        """
        if self._kernel is not None and not recompute_if_exists:
            return self._kernel

        if self._x is None:
            self.get_x()

        print("Computing kernel ...")
        self._kernel = kernel_function(self._x)
        print("Kernel computed")
        return self._kernel

    def save_kernel(self, output_dir):
        """

        Args:
            output_dir:

        Returns:

        """
        if self._kernel is not None:
            filename = path.join(output_dir, 'kernel.txt')
            np.savetxt(filename, self._kernel)
            return filename
        raise Exception("Unable to save the kernel. Kernel must have been computed before.")

    @abc.abstractmethod
    def save_weights_as_nifti(self, weights, output_dir):
        pass


class CAPSVoxelBasedInput(CAPSInput):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, fwhm=0,
                 modulated="on", pvc=None, mask_zeros=True, precomputed_kernel=None):
        """

        Args:
            caps_directory:
            subjects_visits_tsv:
            diagnoses_tsv:
            group_id:
            image_type: 'T1', 'fdg', 'av45', 'pib' or 'flute'
            fwhm:
            modulated:
            mask_zeros:
            precomputed_kernel:
        """

        super(CAPSVoxelBasedInput, self).__init__(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                  image_type, precomputed_kernel=precomputed_kernel)

        self._fwhm = fwhm
        self._modulated = modulated
        self._pvc = pvc
        self._mask_zeros = mask_zeros
        self._orig_shape = None
        self._data_mask = None

        if modulated not in ['on', 'off']:
            raise Exception("Incorrect modulation parameter. It must be one of the values 'on' or 'off'")

    def get_images(self):
        """

        Returns: a list of filenames

        """
        if self._images is not None:
            return self._images

        if self._image_type == 'T1':
            fwhm = '' if self._fwhm == 0 else '_fwhm-%dmm' % int(self._fwhm)

            self._images = [path.join(self._caps_directory, 'subjects', self._subjects[i], self._sessions[i],
                                      't1/spm/dartel/group-' + self._group_id,
                                      '%s_%s_T1w_segm-graymatter_space-Ixi549Space_modulated-%s%s_probability.nii.gz'
                                      % (self._subjects[i], self._sessions[i], self._modulated, fwhm))
                            for i in range(len(self._subjects))]
        else:
            pvc = '' if self._pvc is None else '_pvc-%s' % self._pvc
            fwhm = '' if self._fwhm == 0 else '_fwhm-%dmm' % int(self._fwhm)
            suvr = 'pons' if self._image_type == 'fdg' else 'cerebellumPons'

            self._images = [path.join(self._caps_directory, 'subjects', self._subjects[i], self._sessions[i],
                                      'pet/preprocessing/group-' + self._group_id,
                                      '%s_%s_task-rest_acq-%s_pet_space-Ixi549Space%s_suvr-%s_mask-brain%s_pet.nii.gz'
                                      % (self._subjects[i], self._sessions[i], self._image_type, pvc, suvr, fwhm))
                            for i in range(len(self._subjects))]

        for image in self._images:
            if not path.exists(image):
                raise Exception("File %s doesn't exists." % image)

        return self._images

    def get_x(self):
        """

        Returns: a numpy 2d-array.

        """
        if self._x is not None:
            return self._x

        print('Loading ' + str(len(self.get_images())) + ' subjects')
        self._x, self._orig_shape, self._data_mask = vbio.load_data(self._images, mask=self._mask_zeros)
        print('Subjects loaded')

        return self._x

    def save_weights_as_nifti(self, weights, output_dir):

        if self._images is None:
            self.get_images()

        output_filename = path.join(output_dir, 'weights.nii.gz')
        data = vbio.revert_mask(weights, self._data_mask, self._orig_shape)
        vbio.weights_to_nifti(data, self._images[0], output_filename)


class CAPSRegionBasedInput(CAPSInput):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, atlas,
                 pvc=None, precomputed_kernel=None):
        """

        Args:
            caps_directory:
            subjects_visits_tsv:
            diagnoses_tsv:
            group_id:
            image_type: 'T1', 'fdg', 'av45', 'pib' or 'flute'
            atlas:
            precomputed_kernel:
        """

        super(CAPSRegionBasedInput, self).__init__(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                   image_type, precomputed_kernel=precomputed_kernel)

        self._atlas = atlas
        self._pvc = pvc
        self._orig_shape = None
        self._data_mask = None

        if atlas not in ['AAL2', 'Neuromorphometrics', 'AICHA', 'LPBA40', 'Hammers']:
            raise Exception("Incorrect atlas name. It must be one of the values 'AAL2', 'Neuromorphometrics', 'AICHA', 'LPBA40', 'Hammers' ")

    def get_images(self):
        """

        Returns: a list of filenames

        """
        if self._images is not None:
            return self._images

        if self._image_type == 'T1':
            self._images = [path.join(self._caps_directory, 'subjects', self._subjects[i], self._sessions[i],
                                      't1/spm/dartel/group-' + self._group_id,
                                      'atlas_statistics/', '%s_%s_T1w_space-%s_map-graymatter_statistics.tsv'
                                      % (self._subjects[i], self._sessions[i], self._atlas))
                            for i in range(len(self._subjects))]
        else:
            pvc = '' if self._pvc is None else '_pvc-%s' % self._pvc
            suvr = 'pons' if self._image_type == 'fdg' else 'cerebellumPons'

            self._images = [path.join(self._caps_directory, 'subjects', self._subjects[i], self._sessions[i],
                                      'pet/preprocessing/group-' + self._group_id, 'atlas_statistics',
                                      '%s_%s_task-rest_acq-%s_pet_space-%s%s_suvr-%s_statistics.tsv'
                                      % (self._subjects[i], self._sessions[i], self._image_type, self._atlas, pvc, suvr))
                            for i in range(len(self._subjects))]

        for image in self._images:
            if not path.exists(image):
                raise Exception("File %s doesn't exists." % image)

        return self._images

    def get_x(self):
        """

        Returns: a numpy 2d-array.

        """
        if self._x is not None:
            return self._x

        print('Loading ' + str(len(self.get_images())) + ' subjects')
        self._x = rbio.load_data(self._images, self._subjects)
        print('Subjects loaded')

        return self._x

    def save_weights_as_nifti(self, weights, output_dir):
        """

        Args:
            weights:
            output_dir:

        Returns:

        """

        output_filename = path.join(output_dir, 'weights.nii.gz')
        rbio.weights_to_nifti(weights, self._atlas, output_filename)


class CAPSVertexBasedInput(CAPSInput):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, fwhm, image_type, precomputed_kernel=None):
        super(CAPSVertexBasedInput, self).__init__(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                   image_type, precomputed_kernel)
        self._fwhm = fwhm
        self._image_type = image_type
        self._caps_directory = caps_directory

    def get_images(self):
        import os
        """
        returns list of filnames
        """

        if self._images is not None:
            return self._images

        if self._image_type == 'fdg' and self._images is None:
            self._images = []
            hemi = ['lh', 'rh']
            for i in range(len(self._subjects)):
                self._images.append([os.path.join(self._caps_directory, 'subjects', self._subjects[i], self._sessions[i], 'pet',
                                                  'surface', self._subjects[i] + '_' + self._sessions[i]
                                                  + '_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-' + h
                                                  + '_fwhm-' + str(self._fwhm) + '_projection.mgh') for h in hemi])
            missing_files = []
            missing_files_string_error = ''
            for img in self._images:
                for side in img:
                    if not os.path.exists(side):
                        missing_files.append(side)
                        missing_files_string_error += side + '\n'
            if len(missing_files) > 0:
                raise IOError('Could not find the following files : \n' + missing_files_string_error
                              + '\n' + str(len(missing_files)) + ' files missing')
        return self._images

    def get_x(self):
        from clinica.utils.stream import cprint
        """
        Returns numpy 2D array
        """

        if self._x is not None:
            return self._x

        cprint('Loading ' + str(len(self.get_images())) + ' subjects')
        self._x = vtxbio.load_data(self._images)
        cprint(str(len(self._x)) + ' subjects loaded')
        return self._x

    def save_weights_as_datasurface(self, weights, output_dir):
        import numpy as np
        import nibabel as nib
        import os

        if self._images is None:
            self.get_images()

        sample = nib.load(self._images[0][0])

        infinite_norm = np.max(np.abs(weights))

        left_hemi_data = np.atleast_3d(np.divide(weights[:np.int(weights.size / 2)], infinite_norm))
        left_hemi_mgh = nib.MGHImage(left_hemi_data, affine=sample.affine, header=sample.header)
        nib.save(left_hemi_mgh, os.path.join(output_dir, 'weights_lh.mgh'))

        right_hemi_data = np.atleast_3d(np.divide(weights[np.int(weights.size/2):], infinite_norm))
        right_hemi_mgh = nib.MGHImage(right_hemi_data, affine=sample.affine, header=sample.header)
        nib.save(right_hemi_mgh, os.path.join(output_dir, 'weights_rh.mgh'))
        pass

    def save_weights_as_nifti(self, weights, output_dir):
        pass


class CAPSTSVBasedInput(CAPSInput):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, atlas, dataset,
                 pvc=None, precomputed_kernel=None):
        """

        Args:
            caps_directory:
            subjects_visits_tsv:
            diagnoses_tsv:
            group_id:
            image_type: 'T1', 'fdg', 'av45', 'pib' or 'flute'
            atlas:
            precomputed_kernel:
        """

        super(CAPSTSVBasedInput, self).__init__(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                image_type, precomputed_kernel)

        self._atlas = atlas
        self._pvc = pvc
        self._dataset = dataset

        self._orig_shape = None
        self._data_mask = None

        if atlas not in ['AAL2', 'Neuromorphometrics', 'AICHA', 'LPBA40', 'Hammers']:
            raise Exception("Incorrect atlas name. It must be one of the values 'AAL2', 'Neuromorphometrics', 'AICHA', 'LPBA40', 'Hammers' ")

    def get_images(self):
        """

        Returns: string

        """

        pass

    def get_x(self):
        """

        Returns: a numpy 2d-array.

        """

        # if self._x is not None:
        #    return self._x

        print('Loading TSV subjects')
        string = str('group-' + self._group_id + '_T1w_space-' + self._atlas + '_map-graymatter')

        self._x = tbio.load_data(string, self._caps_directory, self._subjects, self._sessions, self._dataset)

        print('Subjects loaded')

        return self._x

    def save_weights_as_nifti(self, weights, output_dir):
        """

        Args:
            weights:
            output_dir:

        Returns:

        """

        # output_filename = path.join(output_dir, 'weights.nii.gz')

        # rbio.weights_to_nifti(weights, self._atlas, output_filename)
        pass


class CAPSVoxelBasedInputREGSVM(CAPSInput):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, fwhm=0,
                 modulated="on", pvc=None, mask_zeros=True, precomputed_kernel=None):
        """

        Args:
            caps_directory:
            subjects_visits_tsv:
            diagnoses_tsv:
            group_id:
            image_type: 'T1', 'fdg', 'av45', 'pib' or 'flute'
            fwhm:
            modulated:
            mask_zeros:
            precomputed_kernel:
        """

        super(CAPSVoxelBasedInputREGSVM, self).__init__(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                        image_type, precomputed_kernel=precomputed_kernel)

        self._fwhm = fwhm
        self._modulated = modulated
        self._pvc = pvc
        self._mask_zeros = mask_zeros
        self._orig_shape = None
        self._data_mask = None

        if modulated not in ['on', 'off']:
            raise Exception("Incorrect modulation parameter. It must be one of the values 'on' or 'off'")

    def get_images(self):
        """

        Returns: a list of filenames

        """
        if self._images is not None:
            return self._images

        if self._image_type == 'T1':
            fwhm = '' if self._fwhm == 0 else '_fwhm-%dmm' % int(self._fwhm)

            self._images = [path.join(self._caps_directory,
                                      'regul_%s_%s_T1w_segm-graymatter_space-Ixi549Space_modulated-%s%s_probability.nii'
                                      % (self._subjects[i], self._sessions[i], self._modulated, fwhm))
                            for i in range(len(self._subjects))]
        else:
            pvc = '' if self._pvc is None else '_pvc-%s' % self._pvc
            fwhm = '' if self._fwhm == 0 else '_fwhm-%dmm' % int(self._fwhm)
            suvr = 'pons' if self._image_type == 'fdg' else 'cerebellumPons'
            self._images = [path.join(self._caps_directory, 'subjects', self._subjects[i], self._sessions[i],
                                      'pet/preprocessing/group-' + self._group_id,
                                      '%s_%s_task-rest_acq-%s_pet_space-Ixi549Space%s_suvr-%s_mask-brain%s_pet.nii.gz'
                                      % (self._subjects[i], self._sessions[i], self._image_type, pvc, suvr, fwhm))
                            for i in range(len(self._subjects))]

        for image in self._images:
            if not path.exists(image):
                raise Exception("File %s doesn't exists." % image)

        return self._images

    def get_x(self):
        """

        Returns: a numpy 2d-array.

        """
        if self._x is not None:
            return self._x

        print('Loading ' + str(len(self.get_images())) + ' subjects')
        self._x, self._orig_shape, self._data_mask = vbio.load_data(self._images, mask=self._mask_zeros)
        print('Subjects loaded')

        return self._x

    def save_weights_as_nifti(self, weights, output_dir):

        if self._images is None:
            self.get_images()

        output_filename = path.join(output_dir, 'weights.nii.gz')
        data = vbio.revert_mask(weights, self._data_mask, self._orig_shape)
        vbio.weights_to_nifti(data, self._images[0], output_filename)
