# coding: utf8

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "Junhao.Wen@inria.fr"
__status__ = "Development"

from nipype.interfaces import utility as niu
from nipype.interfaces import fsl
from nipype.pipeline import engine as pe
from nipype.workflows.dmri.fsl.utils import (b0_average, compute_readout,)


########################################################################
# NODDI
########################################################################

from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.interfaces.base import (traits, TraitedSpec, File,
                    isdefined)
import os


class EddyNoddiInputSpec(FSLCommandInputSpec):
    in_file = File(exists=True, mandatory=True, argstr='--imain=%s',
                   desc=('File containing all the images to estimate '
                         'distortions for'))
    in_mask = File(exists=True, mandatory=True, argstr='--mask=%s',
                   desc='Mask to indicate brain')
    in_index = File(exists=True, mandatory=True, argstr='--index=%s',
                    desc=('File containing indices for all volumes in --imain '
                          'into --acqp and --topup'))
    in_acqp = File(exists=True, mandatory=True, argstr='--acqp=%s',
                   desc='File containing acquisition parameters')
    in_bvec = File(exists=True, mandatory=True, argstr='--bvecs=%s',
                   desc=('File containing the b-vectors for all volumes in '
                         '--imain'))
    in_bval = File(exists=True, mandatory=True, argstr='--bvals=%s',
                   desc=('File containing the b-values for all volumes in '
                         '--imain'))
    out_base = traits.Str('eddy_corrected', argstr='--out=%s',
                          usedefault=True,
                          desc=('basename for output (warped) image'))
    session = File(exists=True, argstr='--session=%s',
                   desc=('File containing session indices for all volumes in '
                         '--imain'))
    in_topup_fieldcoef = File(exists=True, argstr="--topup=%s",
                              requires=['in_topup_movpar'],
                              desc=('topup file containing the field '
                                    'coefficients'))
    in_topup_movpar = File(exists=True, requires=['in_topup_fieldcoef'],
                           desc='topup movpar.txt file')

    flm = traits.Enum('linear', 'quadratic', 'cubic', argstr='--flm=%s',
                      desc='First level EC model')

    fwhm = traits.Float(desc=('FWHM for conditioning filter when estimating '
                              'the parameters'), argstr='--fwhm=%s')

    niter = traits.Int(5, argstr='--niter=%s', desc='Number of iterations')

    method = traits.Enum('jac', 'lsr', argstr='--resamp=%s',
                         desc=('Final resampling method (jacobian/least '
                               'squares)'))
    repol = traits.Bool(False, argstr='--repol',
                        desc='Detect and replace outlier slices')
    num_threads = traits.Int(1, usedefault=True, nohash=True,
                             desc="Number of openmp threads to use")
    data_is_shelled = traits.Bool(False, argstr='--data_is_shelled',
                        desc='Bypass any checking and eddy will proceed as if data was shelled')


class EddyNoddiOutputSpec(TraitedSpec):
    out_corrected = File(exists=True,
                         desc=('4D image file containing all the corrected '
                               'volumes'))
    out_parameter = File(exists=True,
                         desc=('text file with parameters definining the '
                               'field and movement for each scan'))


class EddyNoddi(FSLCommand):
    """
    This is a modified version of nipype EDDY interface which targets at the preprocessing of NODDI data(opposite two phase encoding direction dwis)
    Interface for FSL eddy, a tool for estimating and correcting eddy
    currents induced distortions. `User guide
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Eddy/UsersGuide>`_ and
    `more info regarding acqp file
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_do_I_know_what_to_put_into_my_--acqp_file>`_.

    Note:
        If you run the versions before 5.0.10, please know that we use another version of binary of eddy in FSL and also, we set these two flags to be true:
            ecc.inputs.repol = True
            ecc.inputs.data_is_shelled = True
            Also, need to install eddy_cuda7.5 on your linux machine, and eddy on Mac.
        If your FSL version is 5.0.10, everything is OK


    Examples
    --------

    from nipype.interfaces.fsl import Eddy
    eddy = Eddy()
    eddy.inputs.in_file = 'epi.nii'
    eddy.inputs.in_mask  = 'epi_mask.nii'
    eddy.inputs.in_index = 'epi_index.txt'
    eddy.inputs.in_acqp  = 'epi_acqp.txt'
    eddy.inputs.in_bvec  = 'bvecs.scheme'
    eddy.inputs.in_bval  = 'bvals.scheme'
    eddy.cmdline #doctest: +ELLIPSIS
    'eddy --acqp=epi_acqp.txt --bvals=bvals.scheme --bvecs=bvecs.scheme \
--imain=epi.nii --index=epi_index.txt --mask=epi_mask.nii \
--out=.../eddy_corrected'
    res = eddy.run() # doctest: +SKIP

    """
    _cmd = 'eddy'
    input_spec = EddyNoddiInputSpec
    output_spec = EddyNoddiOutputSpec

    _num_threads = 1

    def __init__(self, **inputs):
        super(EddyNoddi, self).__init__(**inputs)
        self.inputs.on_trait_change(self._num_threads_update, 'num_threads')

        if not isdefined(self.inputs.num_threads):
            self.inputs.num_threads = self._num_threads
        else:
            self._num_threads_update()

    def _num_threads_update(self):
        self._num_threads = self.inputs.num_threads
        if not isdefined(self.inputs.num_threads):
            if 'OMP_NUM_THREADS' in self.inputs.environ:
                del self.inputs.environ['OMP_NUM_THREADS']
        else:
            self.inputs.environ['OMP_NUM_THREADS'] = str(self.inputs.num_threads)

    def _format_arg(self, name, spec, value):
        if name == 'in_topup_fieldcoef':
            return spec.argstr % value.split('_fieldcoef')[0]
        if name == 'out_base':
            return spec.argstr % os.path.abspath(value)
        return super(EddyNoddi, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_corrected'] = os.path.abspath('%s.nii.gz' % self.inputs.out_base)
        outputs['out_parameter'] = os.path.abspath('%s.eddy_parameters' % self.inputs.out_base)
        return outputs


def sdc_peb_noddi(name='sdc_ped_noddi',
            epi_params=dict(echospacing=0.77e-3,
                            acc_factor=3,
                            enc_dir='y-',
                            epi_factor=1),
            alt_epi_params=dict(echospacing=0.77e-3,
                               acc_factor=3,
                               enc_dir='y',
                               epi_factor=1)):
    """
    SDC stands for susceptibility distortion correction. PEB stands for
    phase-encoding-based.

    The phase-encoding-based (PEB) method implements SDC by acquiring
    diffusion images with two different enconding directions [Andersson2003]_.
    The most typical case is acquiring with opposed phase-gradient blips
    (e.g. *AP* and *PA*, or equivalently, *-y* and *y*)
    as in [Chiou2000]_, but it is also possible to use orthogonal
    configurations [Cordes2000]_ (e.g. *AP* and *LR*,
    or equivalently *-y* and *x*).
    This workflow uses the implementation of FSL
    (`TOPUP <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP>`_).

    Example
    -------

     from nipype.workflows.dmri.fsl.artifacts import sdc_peb

    .. admonition:: References

      .. [Andersson2003] Andersson JL et al., `How to correct susceptibility
        distortions in spin-echo echo-planar images: application to diffusion
        tensor imaging <http://dx.doi.org/10.1016/S1053-8119(03)00336-7>`_.
        Neuroimage. 2003 Oct;20(2):870-88. doi: 10.1016/S1053-8119(03)00336-7

      .. [Cordes2000] Cordes D et al., Geometric distortion correction in EPI
        using two images with orthogonal phase-encoding directions, in Proc.
        ISMRM (8), p.1712, Denver, US, 2000.

      .. [Chiou2000] Chiou JY, and Nalcioglu O, A simple method to correct
        off-resonance related distortion in echo planar imaging, in Proc.
        ISMRM (8), p.1712, Denver, US, 2000.

    """

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_b0s', 'in_full_file', 'in_full_bvals']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_file', 'out_vsm', 'out_movpar', 'out_field_hz', 'out_enc_file', 'out_mask']), name='outputnode')

    avg_b0 = pe.Node(niu.Function(
        input_names=['in_dwi', 'in_bval'], output_names=['out_file'],
        function=b0_average), name='b0_avg')

    topup = pe.Node(fsl.TOPUP(), name='topup')
    # topup.inputs.encoding_direction = [epi_params['enc_dir'],
    #                                    alt_epi_params['enc_dir']]

    readout = compute_readout(epi_params)
    readout_alt = compute_readout(alt_epi_params)
    # topup.inputs.readout_times = [readout,
    #                               readout_alt]

    topup_acq = pe.Node(niu.Function(input_names=['in_file', 'epi_params', 'alt_epi_params', 'readout', 'readout_alt'],
                                     output_names=['out_file'], function=gen_acq_noddi), name='generate_acq_txt_topup')
    topup_acq.inputs.epi_params = epi_params
    topup_acq.inputs.alt_epi_params = alt_epi_params
    topup_acq.inputs.readout = readout
    topup_acq.inputs.readout_alt = readout_alt

    undis_mask = pe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                       name='mask_from_topup')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, topup, [('in_b0s', 'in_file')]),
        (inputnode, avg_b0, [('in_full_file', 'in_dwi'),
                             ('in_full_bvals', 'in_bval')]),
        (avg_b0, undis_mask, [('out_file', 'in_file')]),
        (inputnode, topup_acq, [('in_b0s', 'in_file')]),
        (topup_acq, topup, [('out_file', 'encoding_file')]),
        (topup, outputnode, [('out_fieldcoef', 'out_vsm')]),
        (topup, outputnode, [('out_movpar', 'out_movpar')]),
        (topup, outputnode, [('out_corrected', 'out_file')]),
        (topup, outputnode, [('out_field', 'out_field_hz')]),
        (topup_acq, outputnode, [('out_file', 'out_enc_file')]),
        (undis_mask, outputnode, [('mask_file', 'out_mask')])
    ])
    return wf

################################################################################
# Utilities for the pipeline
################################################################################


def gen_index_noddi(in_bval, b0_index):
    """
    This is a function to generate the index file for FSL eddy
    :param in_bval:
    :param b0_index:
    :return:
    """
    import numpy as np
    import os
    out_file = os.path.abspath('index.txt')
    bvals = np.loadtxt(in_bval)
    vols = len(bvals)
    index_list = []
    for i in range(0, len(b0_index)):
        if i == (len(b0_index) - 1):
            index_list.extend([i+1] * (vols - b0_index[i]))
        else:
            index_list.extend([i+1] * (b0_index[i+1] - b0_index[i]))
    index_array = np.asarray(index_list)
    try:
        len(index_list) == vols
    except ValueError:
        print "It seems that you do not define the index file for FSL eddy correctly!"
    np.savetxt(out_file, index_array.T)

    return out_file


def gen_acq_noddi(in_file, epi_params, alt_epi_params, readout, readout_alt):
    """
    This is a function to generate the FSL topup acq.txt file
    :param in_file:
    :param epi_params:
    :param alt_epi_params:
    :param readout:
    :param readout_alt:
    :return:
    """
    import numpy as np
    import os
    import nibabel as nb
    out_file = os.path.abspath('acq.txt')
    vols = nb.load(in_file).get_data().shape[-1]
    arr = np.ones([vols, 4])
    for i in range(vols):
        if i < vols/2:
            if epi_params['enc_dir'] == 'y-':
                arr[i, :] = np.array((0, -1, 0, readout))
            elif epi_params['enc_dir'] == 'y':
                arr[i, :] = np.array((0, 1, 0, readout))
            elif epi_params['enc_dir'] == 'x':
                arr[i, :] = np.array((0, 1, 0, readout))
            elif epi_params['enc_dir'] == 'x-':
                arr[i, :] = np.array((0, -1, 0, readout))
            elif epi_params['enc_dir'] == 'z':
                arr[i, :] = np.array((0, 1, 0, readout))
            elif epi_params['enc_dir'] == 'x-':
                arr[i, :] = np.array((0, -1, 0, readout))
        else:
            if alt_epi_params['enc_dir_alt'] == 'y-':
                arr[i, :] = np.array((0, -1, 0, readout_alt))
            elif alt_epi_params['enc_dir_alt'] == 'y':
                arr[i, :] = np.array((0, 1, 0, readout_alt))
            elif alt_epi_params['enc_dir_alt'] == 'x':
                arr[i, :] = np.array((0, 1, 0, readout_alt))
            elif alt_epi_params['enc_dir_alt'] == 'x-':
                arr[i, :] = np.array((0, -1, 0, readout_alt))
            elif alt_epi_params['enc_dir_alt'] == 'z':
                arr[i, :] = np.array((0, 1, 0, readout_alt))
            elif alt_epi_params['enc_dir_alt'] == 'x-':
                arr[i, :] = np.array((0, -1, 0, readout_alt))
    np.savetxt(out_file, arr)

    return out_file


def merge_noddi_ped(in_file, in_bvec, in_bval, alt_file, alt_bvec, alt_bval):
    """
    This is to merge the two ped images and also concatenate the bvecs and bvals
    :return:
    """
    from clinica.utils.dwi import merge_volumes_tdim
    import os.path as op
    import os

    out_bvals_tab = op.abspath('merged.bval')
    out_bvals = op.abspath('merged_bvals.bval')
    cmd_bval = 'paste ' + in_bval + ' ' + alt_bval + ' > ' + out_bvals_tab
    os.system(cmd_bval)
    with open(out_bvals_tab, 'r+') as fin, open(out_bvals, 'w') as fout:
        for line in fin:
            fout.write(line.replace('\t', ' '))

    out_bvecs_tab = op.abspath('merged.bvec')
    out_bvecs = op.abspath('merged_bvecs.bvec')
    cmd_bvec = 'paste ' + in_bvec + ' ' + alt_bvec + ' > ' + out_bvecs_tab
    os.system(cmd_bvec)
    with open(out_bvecs_tab, 'r+') as fin, open(out_bvecs, 'w') as fout:
        for line in fin:
            fout.write(line.replace('\t', ' '))
    out_file = merge_volumes_tdim(in_file, alt_file)

    return out_file, out_bvals, out_bvecs


def grab_noddi_bids_files(bids_directory, tsv):
    """
    This is a function to get all the files that need for noddi preprocessing in BIDS structure
    :param bids_directory:
    :param tsv:
    :return:
    """
    import os, csv

    subject_list = []
    session_list = []
    with open(tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            if row[0] == 'participant_id':
                continue
            else:
                subject_list.append(row[0])
                session_list.append(row[1])
        bids_directory = os.path.expanduser(bids_directory)  # change the relative path to be absolute path

    bids_ap_dwi = []
    bids_ap_dwi_bvec = []
    bids_ap_dwi_bval = []

    bids_pa_dwi = []
    bids_pa_dwi_bvec = []
    bids_pa_dwi_bval = []

    # the number of subject_list and session_list should be the same
    try:
        len(subject_list) == len(session_list)
    except RuntimeError:
        print "It seems that the nuber of session_list and subject_list are not in the same length, please check"
        raise

    num_subject = len(subject_list)
    for i in xrange(num_subject):
        # AP
        subjec_nii = os.path.join(bids_directory, subject_list[i], session_list[i], 'dwi', subject_list[i] + '_' + session_list[i] + '_seq-mshellAP_dwi.nii.gz')
        bids_ap_dwi += [subjec_nii]

        subjec_bvec = os.path.join(bids_directory, subject_list[i], session_list[i], 'dwi', subject_list[i] + '_' + session_list[i] + '_seq-mshellAP_dwi.bvec')
        bids_ap_dwi_bvec += [subjec_bvec]

        subjec_bval = os.path.join(bids_directory, subject_list[i], session_list[i], 'dwi', subject_list[i] + '_' + session_list[i] + '_seq-mshellAP_dwi.bval')
        bids_ap_dwi_bval += [subjec_bval]

        # PA
        subjec_nii = os.path.join(bids_directory, subject_list[i], session_list[i], 'dwi',
                                  subject_list[i] + '_' + session_list[i] + '_seq-mshellPA_dwi.nii.gz')
        bids_pa_dwi += [subjec_nii]

        subjec_bvec = os.path.join(bids_directory, subject_list[i], session_list[i], 'dwi',
                                   subject_list[i] + '_' + session_list[i] + '_seq-mshellPA_dwi.bvec')
        bids_pa_dwi_bvec += [subjec_bvec]

        subjec_bval = os.path.join(bids_directory, subject_list[i], session_list[i], 'dwi',
                                   subject_list[i] + '_' + session_list[i] + '_seq-mshellPA_dwi.bval')
        bids_pa_dwi_bval += [subjec_bval]

    return bids_ap_dwi, bids_ap_dwi_bvec, bids_ap_dwi_bval, bids_pa_dwi, bids_pa_dwi_bvec, bids_pa_dwi_bval


def get_subid_sesid(in_file, caps_directory):
    """
    This is to extract the base_directory for the DataSink including participant_id and sesion_id in CAPS directory, also the tuple_list for substitution
    :param in_file:
    :return: base_directory for DataSink
    """
    import os

    identifier = (in_file.split('/')[-1]).split('_seq-mshell')[0]
    participant_id = identifier.split('_')[0]
    session_id = identifier.split('_')[1]
    base_directory = os.path.join(caps_directory, 'subjects', participant_id, session_id, 'dwi', 'preprocessing')
    subst_tuple_list = [
        ('eddy_corrected_avg_b0_brain_mask.nii.gz', identifier + '_dwi_space-b0_brainmask.nii.gz'),
        ('merged_bvecs_rotated.bvec', identifier + '_dwi_space-b0_preproc.bvec'),
        ('merged_bvals.bval', identifier + '_dwi_space-b0_preproc.bval'),
        ('merged_files.nii.gz', identifier + '_dwi_space-b0_preproc.nii.gz'),

        # ('acq.txt', identifier + '_acq.txt'),
        # ('b0_base_fieldcoef.nii.gz', identifier + '_fieldcoef.nii.gz'),
        # ('b0_base_movpar.txt', identifier + '_movpar.txt'),
        # ('b0_corrected.nii.gz', identifier + '_corrected_unwarped.nii.gz'),
        # ('b0_field.nii.gz', identifier + '_field_hz.nii.gz'),
        # ('eddy_corrected.eddy_parameters', identifier + '_eddy_corrected.eddy_parameters'),
        # ('eddy_corrected.nii.gz', identifier + '_eddy_corrected.nii.gz')
    ]

    return base_directory, subst_tuple_list


def noddi_preprocessing_twoped(caps_directory, name='noddi_preprocessing_topup_eddy',
                     epi_params=dict(echo_spacing=0.77e-3,
                                     acc_factor=3,
                                     enc_dir='y-',
                                     epi_factor=1),
                     alt_epi_params=dict(echo_spacing=0.77e-3,
                                        acc_factor=3,
                                        enc_dir='y',)):
    """
    This is a preprocessing pipeline for typical NODDI data with two phase encoding direction(ap/j-/y-, pa/j/y) images which will implement the new tool
    from FSL (topup and eddy) to correct the susceptibility-induced off-resonance field distortion by topup, eddy current-induced off-resonance
    field correction and also head motion correction between volume by EDDY(new version of eddy_correct).
    WARNING: Make sure that for your analysis, the protocol of all the subjects(parameters) should be the same.

    The steps are:
    1) data preparation, extract the B0 and concatenate them in one single 4D file, and concatenate the two opposite ped images to one singe dwi(ap+pa),
        also for the bvec and bval files
    2) run topup just for the concatenated B0 images.
    3) run eddy with the previous topup output for the concatenated dwi,
        note that you need to run eddy with the option  --data_is_shelled   (you need a different binary than the current fsl release),
        download: https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/ ;
        explanation: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;d3bf61ea.1610
    :param caps_directory:
    :param name:
    :param epi_params:
    :param alt_epi_params:
    :return:

    Note: for PREVDEMALS ICM, the data:
        DTI_1_ap d66_b2200 (60 dwis + 6 b0)
        DTI_3_ap d38_b700  (32 dwis + 6 b0)
        DTI_5_ap d10_b300  (9 dwis + 1 b0)
        DTI_2_pa d66_b2200
        DTI_4_pa d38_b700
        DTI_6_pa d10_b300

        PREVDEMALS ICM has two json files in BIDS for each subject which contains all the parameters for this sequence, please check out the parameters there.

    """

    from clinica.utils.dwi import b0_dwi_split
    # import merge_noddi_ped, gen_index_noddi, EddyNoddi, sdc_peb_noddi, get_subid_sesid
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe
    from nipype.workflows.dmri.fsl.utils import (b0_indices, b0_average, eddy_rotate_bvecs)
    from nipype.interfaces import fsl
    import nipype.interfaces.io as nio

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_file', 'in_bvec', 'in_bval', 'alt_file', 'alt_bvec', 'alt_bval']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['ecc_out_file', 'out_mask', 'out_bvec', 'original_merged_bval']), name='outputnode')

    merge_two_ped = pe.Node(niu.Function(input_names=['in_file', 'in_bvec', 'in_bval', 'alt_file', 'alt_bvec', 'alt_bval'],
                                         output_names=['out_file', 'out_bvals', 'out_bvecs'], function=merge_noddi_ped), name='merge_ped_images')

    extract_b0_dwis = pe.Node(
            niu.Function(
                input_names=['in_dwi', 'in_bval', 'in_bvec'],
                output_names=['out_b0', 'out_dwi', 'out_bvals', 'out_bvecs'],
                function=b0_dwi_split),
            name='extract_dwis_b0')

    list_b0 = pe.Node(niu.Function(
        input_names=['in_bval'], output_names=['out_idx'],
        function=b0_indices), name='find_b0_indices')

    generate_index_eddy = pe.Node(niu.Function(input_names=['in_bval', 'b0_index'],
        output_names=['eddy_index'],
        function=gen_index_noddi), name='generate_index_eddy')

    avg_b0 = pe.Node(niu.Function(
        input_names=['in_dwi', 'in_bval'], output_names=['out_file'],
        function=b0_average), name='b0_avg_post')

    bet_dwi = pe.Node(
            fsl.BET(frac=0.5, mask=True, robust=True),
            name='bet_dwi_post')

    sdc = sdc_peb_noddi(epi_params=epi_params, alt_epi_params=alt_epi_params)

    ecc = pe.Node(EddyNoddi(method='jac'), name='fsl_eddy_noddi')
    ecc.inputs.repol = True
    ecc.inputs.data_is_shelled = True

    rot_bvec = pe.Node(niu.Function(
        input_names=['in_bvec', 'eddy_params'], output_names=['out_file'],
        function=eddy_rotate_bvecs), name='rotate_bvec')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, merge_two_ped, [('in_file', 'in_file'),
                               ('in_bval', 'in_bval'),
                               ('in_bvec', 'in_bvec'),
                               ('alt_file', 'alt_file'),
                               ('alt_bvec', 'alt_bvec'),
                               ('alt_bval', 'alt_bval')]),
        (merge_two_ped, extract_b0_dwis,  [('out_file', 'in_dwi'),
                                           ('out_bvals', 'in_bval'),
                                           ('out_bvecs', 'in_bvec')]),
        (extract_b0_dwis, sdc,  [('out_b0', 'inputnode.in_b0s')]),
        (merge_two_ped, sdc,  [('out_file', 'inputnode.in_full_file')]),
        (merge_two_ped, sdc,  [('out_bvals', 'inputnode.in_full_bvals')]),
        (sdc, ecc,    [('outputnode.out_mask', 'in_mask'),
                       ('outputnode.out_enc_file', 'in_acqp'),
                       ('topup.out_fieldcoef', 'in_topup_fieldcoef'),
                       ('topup.out_movpar', 'in_topup_movpar')]),
        (merge_two_ped, list_b0, [('out_bvals', 'in_bval')]),
        (list_b0, generate_index_eddy, [('out_idx', 'b0_index')]),
        (merge_two_ped, generate_index_eddy, [('out_bvals', 'in_bval')]),
        (generate_index_eddy, ecc, [('eddy_index', 'in_index')]),
        (merge_two_ped, ecc, [('out_file', 'in_file'),
                              ('out_bvals', 'in_bval'),
                              ('out_bvecs', 'in_bvec')]),
        (merge_two_ped, rot_bvec, [('out_bvecs', 'in_bvec')]),
        (ecc, rot_bvec, [('out_parameter', 'eddy_params')]),
        (ecc,       avg_b0,   [('out_corrected', 'in_dwi')]),
        (inputnode, avg_b0,   [('in_bval', 'in_bval')]),
        (avg_b0,  bet_dwi,   [('out_file', 'in_file')]),
        # orginal files
        (merge_two_ped,    outputnode, [('out_bvals', 'original_merged_bval')]),
        # ecc files
        (ecc, outputnode, [('out_corrected', 'ecc_out_file')]),
        (rot_bvec, outputnode, [('out_file', 'out_bvec')]),
        (bet_dwi, outputnode, [('mask_file', 'out_mask')])
    ])
    return wf


def epi_params(echo_spacing, acc_factor, enc_dir, enc_dir_alt, epi_factor):
    """
    This function is to create a dict based on the input for epi parameters for both phase encoding directions
    :param echo_spacing:
    :param acc_factor:
    :param enc_dir:
    :param enc_dir_alt:
    :param epi_factor:
    :return:
    """
    epi_param, epi_param_alt = dict(echospacing=echo_spacing, acc_factor=acc_factor, enc_dir=enc_dir, epi_factor=epi_factor), dict(echospacing=echo_spacing, acc_factor=acc_factor, enc_dir=enc_dir_alt, epi_factor=epi_factor)
    return epi_param, epi_param_alt


def create_list_tuple(bids_ap_dwi, bids_ap_dwi_bvec, bids_ap_dwi_bval, bids_pa_dwi, bids_pa_dwi_bvec, bids_pa_dwi_bval):
    list_tuple = []
    for i in xrange(len(bids_ap_dwi)):
        ls = []
        ls.append(bids_ap_dwi[i])
        ls.append(bids_ap_dwi_bvec[i])
        ls.append(bids_ap_dwi_bval[i])
        ls.append(bids_pa_dwi[i])
        ls.append(bids_pa_dwi_bvec[i])
        ls.append(bids_pa_dwi_bval[i])
        ls = tuple(ls)
        list_tuple.append(ls)
    return list_tuple

