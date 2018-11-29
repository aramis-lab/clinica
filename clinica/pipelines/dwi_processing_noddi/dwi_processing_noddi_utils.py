# coding: utf8


def amico_noddi(root_dir, subject_id, in_bval, in_bvec, noddi_preprocessed, mask_file, bStep):
    """
     This is a function to call AMICO based on amico dataset structure

    Args:
        root_dir: path to the folder containing all the subjects from one study
        subject_id: path to the subject id
        in_bval: path to the bval
        in_bvec: paht to the bvec
        noddi_preprocessed: path to the preprocessed NODDI mri
        mask_file: path to the mask file
        bStep:

    Returns:

    """

    import amico, os
    from os.path import join as opj
    from os import chdir as oc

    # checking if the files exist
    try:
        os.path.isfile(in_bval) and os.path.isfile(in_bvec) and os.path.isfile(noddi_preprocessed) and os.path.isfile(mask_file)
    except OSError:
        raise

    # Get the name of the bvec and bval files
    bval_name = in_bval.split('/')[-1]
    bvec_name = in_bvec.split('/')[-1]
    noddi_preprocessed = noddi_preprocessed.split('/')[-1]
    mask_file = mask_file.split('/')[-1]
    # Initial the amico setup config
    amico.core.setup()

    ae = amico.Evaluation(root_dir, subject_id)
    subject_path = opj(root_dir, subject_id)

    scheme_filename = amico.util.fsl2scheme(opj(subject_path, bval_name), opj(subject_path, bvec_name), bStep=bStep)
    oc(subject_path)
    # load the data
    ae.load_data(dwi_filename=noddi_preprocessed, scheme_filename=scheme_filename, mask_filename=mask_file)

    # Compute the response functions
    ae.set_model("NODDI")
    ae.generate_kernels()
    ae.load_kernels()

    # fit the model
    ae.fit()
    ae.save_results()

    # get the output file for datasinker
    FIT_ICVF = os.path.join(subject_path, 'AMICO', 'NODDI', 'FIT_ICVF.nii.gz')
    FIT_ISOVF = os.path.join(subject_path, 'AMICO', 'NODDI', 'FIT_ISOVF.nii.gz')
    FIT_OD = os.path.join(subject_path, 'AMICO', 'NODDI', 'FIT_OD.nii.gz')

    return FIT_ICVF, FIT_ISOVF, FIT_OD, root_dir, subject_id


def create_amico_structure(working_directory, subject_id, noddi_preprocessed, in_bval, in_bvec, mask_file):
    """
        This is to create a AMICO data structure from the original CAPS structure for the preprocessed data which will be ready to fit AMICO noddi model
    Args:
        working_directory:
        subject_id:
        noddi_preprocessed:
        in_bval:
        in_bvec:
        mask_file:

    Returns:

    """

    import os
    import errno
    from shutil import copyfile as cp

    root_dir = working_directory

    try:
        os.makedirs(os.path.join(root_dir, subject_id))
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    cp(noddi_preprocessed, os.path.join(root_dir, subject_id, 'Noddi.nii.gz'))
    cp(in_bval, os.path.join(root_dir, subject_id, 'Noddi.bval'))
    cp(in_bvec, os.path.join(root_dir, subject_id, 'Noddi.bvec'))
    cp(mask_file, os.path.join(root_dir, subject_id, 'Noddi_mask.nii.gz'))

    noddi_preprocessed = os.path.join(root_dir, subject_id, 'Noddi.nii.gz')
    in_bval = os.path.join(root_dir, subject_id, 'Noddi.bval')
    in_bvec = os.path.join(root_dir, subject_id, 'Noddi.bvec')
    mask_file = os.path.join(root_dir, subject_id, 'Noddi_mask.nii.gz')

    return root_dir, subject_id, noddi_preprocessed, in_bval, in_bvec, mask_file



def grab_noddi_preprocessed_files(caps_directory, tsv):
    """
        This is a function to grab the files from the preprocessed data
    Args:
        caps_directory:
        tsv:

    Returns:

    """
    import os, csv
    import tempfile
    from time import time, strftime, localtime
    import clinica.iotools.utils.data_handling as cdh

    if tsv == None:
        output_dir = tempfile.mkdtemp()
        timestamp = strftime('%Y%m%d_%H%M%S', localtime(time()))
        tsv_file = '%s_subjects_sessions_list.tsv' % timestamp
        tsv = os.path.join(output_dir, tsv_file)

        cdh.create_subs_sess_list(
            input_dir=caps_directory,
            output_dir=output_dir,
            file_name=tsv_file,
            is_bids_dir=False)

    subject_list = []
    session_list = []
    subject_id_list = []
    with open(tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            if row[0] == 'participant_id':
                continue
            else:
                subject_list.append(row[0])
                session_list.append(row[1])
                subject_id_list.append((row[0]+'_'+row[1]))
        caps_directory = os.path.expanduser(caps_directory)


    noddi_preprocessed_dwi = []
    noddi_preprocessed_bvec = []
    noddi_preprocessed_bval = []
    noddi_preprocessed_mask = []


    # the number of subject_list and session_list should be the same
    try:
        len(subject_list) == len(session_list)
    except RuntimeError:
        print "It seems that the nuber of session_list and subject_list are not in the same length, please check"
        raise

    num_subject = len(subject_list)
    for i in xrange(num_subject):
        # AP
        subject_nii = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_dwi_space-b0_preproc.nii.gz')
        noddi_preprocessed_dwi += [subject_nii]

        subject_bvec = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_dwi_space-b0_preproc.bvec')
        noddi_preprocessed_bvec += [subject_bvec]

        subject_bval = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_dwi_space-b0_preproc.bval')
        noddi_preprocessed_bval += [subject_bval]

        subject_mask = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_dwi_space-b0_brainmask.nii.gz')
        noddi_preprocessed_mask += [subject_mask]

    return subject_id_list, noddi_preprocessed_dwi, noddi_preprocessed_bvec, noddi_preprocessed_bval, noddi_preprocessed_mask


def grab_noddi_preprocessed_files_oneshell_adni(caps_directory, tsv):
    """
        This is a function to grab the files from the preprocessed data
    Args:
        caps_directory:
        tsv:

    Returns:

    """
    import os, csv

    subject_list = []
    session_list = []
    subject_id_list = []
    with open(tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            if row[0] == 'participant_id':
                continue
            else:
                subject_list.append(row[0])
                session_list.append(row[1])
                subject_id_list.append((row[0]+'_'+row[1]))
        caps_directory = os.path.expanduser(caps_directory)


    noddi_preprocessed_dwi = []
    noddi_preprocessed_bvec = []
    noddi_preprocessed_bval = []
    noddi_preprocessed_mask = []


    # the number of subject_list and session_list should be the same
    try:
        len(subject_list) == len(session_list)
    except RuntimeError:
        print "It seems that the nuber of session_list and subject_list are not in the same length, please check"
        raise

    num_subject = len(subject_list)
    for i in xrange(num_subject):
        # AP
        subject_nii = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_acq-axial_dwi_space-T1w_preproc.nii.gz')
        noddi_preprocessed_dwi += [subject_nii]

        subject_bvec = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing',  subject_id_list[i] + '_acq-axial_dwi_space-T1w_preproc.bvec')
        noddi_preprocessed_bvec += [subject_bvec]

        subject_bval = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_acq-axial_dwi_space-T1w_preproc.bval')
        noddi_preprocessed_bval += [subject_bval]

        subject_mask = os.path.join(caps_directory, 'subjects', subject_list[i], session_list[i], 'dwi', 'preprocessing', subject_id_list[i] + '_acq-axial_dwi_space-T1w_brainmask.nii.gz')
        noddi_preprocessed_mask += [subject_mask]

    return subject_id_list, noddi_preprocessed_dwi, noddi_preprocessed_bvec, noddi_preprocessed_bval, noddi_preprocessed_mask


def delete_amico(root_dir, subject_id):
    """
        This is a function to delete the temporal folder containing Amico structure after running the pipeline to save space,
        but this will just delete the amico structure images, not the noddi maps results
    Args:
        root_dir:

    Returns:

    """
    import os
    from glob import glob

    files_list = glob(os.path.join(os.path.expanduser(root_dir), subject_id, 'Noddi*'))
    if len(files_list) == 0:
        pass
    else:
        [os.remove(f) for f in files_list]


def runmatlab(output_dir, noddi_img, brain_mask, roi_mask, bval, bvec, prefix, bStep, num_cores, path_to_matscript, noddi_toolbox_dir, nifti_matlib_dir):
    """
    The wrapper to call noddi matlab script.
    Args:
        output_dir:
        noddi_img:
        brain_mask:
        roi_mask:
        bval:
        bvec:
        prefix:
        bStep:
        num_cores:

    Returns:

    """
    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
    from os.path import join
    import sys, os
    # here, we check out the os, basically, clinica works for linux and MAC OS X.
    if sys.platform.startswith('linux'):
        print "###Note: your platform is linux, the default command line for Matlab(matlab_cmd) is matlab, but you can also export a variable MATLABCMD,  which points to your matlab,  in your .bashrc to set matlab_cmd, this can help you to choose which Matlab to run when you have more than one Matlab. "
    elif sys.platform.startswith('darwin'):
        try:
            if not 'MATLABCMD' in os.environ:
                raise RuntimeError(
                    "###Note: your platform is MAC OS X, the default command line for Matlab(matlab_cmd) is matlab, but it does not work on OS X, you mush export a variable MATLABCMD, which points to your matlab, in your .bashrc to set matlab_cmd. Note, Mac os x will always choose to use OpengGl hardware mode.")
        except Exception as e:
            print(str(e))
            exit(1)
    else:
        print "Clinica will not work on your platform "

    MatlabCommand.set_default_matlab_cmd(
        get_matlab_command())  # this is to set the matlab_path(os.environ) in your bashrc file, to choose which version of matlab do you wanna use
    # here, set_default_matlab_cmd is a @classmethod
    matlab = MatlabCommand()

    # add the dynamic traits
    # openGL_trait = traits.Bool(True, argstr='-nosoftwareopengl', usedefault=True, desc='Switch on hardware openGL', nohash=True)
    # matlab.input_spec.add_trait(matlab.input_spec(), 'nosoftwareopengl', openGL_trait() )
    if sys.platform.startswith('linux'):
        matlab.inputs.args = '-nosoftwareopengl'  # Bug, for my laptop, it does not work, but the command line does have the flag -nosoftwareopengl, we should try on other computer's matlab to check if this flag works!
    matlab.inputs.paths = path_to_matscript  # CLINICA_HOME, this is the path to add into matlab, addpath

    matlab.inputs.script = """
    noddiprocessing('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %d');
    """ % (output_dir, noddi_img, brain_mask, roi_mask, bval, bvec, prefix, bStep, noddi_toolbox_dir, nifti_matlib_dir, num_cores)  # here, we should define the inputs for the matlab function that you want to use
    matlab.inputs.mfile = True  # this will create a file: pyscript.m , the pyscript.m is the default name
    matlab.inputs.single_comp_thread = False  # this will stop runing with single thread
    matlab.inputs.logfile = join(output_dir, prefix + "_matlab_output.log")
    print "Matlab logfile is located in the folder: %s" % matlab.inputs.logfile
    print "Matlab script command = %s" % matlab.inputs.script
    print "MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread
    print "MatlabCommand choose which matlab to use(matlab_cmd): %s" % get_matlab_command()
    if sys.platform.startswith('linux'):
        print "MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args

    matlab.run()

    # grab the output images
    fit_icvf = os.path.join(output_dir, prefix+'_ficvf.nii')
    fit_isovf = os.path.join(output_dir, prefix+'_fiso.nii')
    fit_od = os.path.join(output_dir, prefix+'_odi.nii')

    return fit_icvf, fit_isovf, fit_od


def get_subid_sesid(in_file, caps_directory):
    """
    This is to extract the base_directory for the DataSink including participant_id and sesion_id in CAPS directory, also the tuple_list for substitution
    :param subject_id:
    :return: base_directory for DataSink
    """
    import os

    identifier = (in_file.split('/')[-1]).split('_ficvf')[0]
    participant_id = identifier.split('_')[0]
    session_id = identifier.split('_')[1]
    base_directory = os.path.join(caps_directory, 'subjects', participant_id, session_id, 'dwi', 'noddi_based_processing', 'native_space')

    subst_tuple_list = [
        (identifier + '_ficvf.nii.gz', identifier + '_space-b0_NDI.nii.gz'),
        (identifier + '_fiso.nii.gz', identifier + '_space-b0_FWF.nii.gz'),
        (identifier + '_odi.nii.gz', identifier + '_space-b0_ODI.nii.gz'),
        (r'/trait_added/', r''),
        (r'/processing_datasinker\d{1,4}/', r''),
        (r'/fit_icvf\d{1,4}/', r''),
        (r'/fit_isovf\d{1,4}/', r''),
        (r'/fit_od\d{1,4}/', r''),
    ]

    return base_directory, subst_tuple_list


def compress_nii(in_file, same_dir=True):
    """
        This is a function to compress the resulting nii images
    Args:
        in_file:
        same_dir:

    Returns:

    """
    from os import getcwd, remove
    from os.path import abspath, join
    import gzip
    import shutil
    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, basestring):  # type(in_file) is list:
        return [compress_nii(f, same_dir) for f in in_file]

    orig_dir, base, ext = split_filename(str(in_file))

    # Already compressed
    if ext[-3:].lower() == ".gz":
        return in_file

    # Not compressed
    if same_dir:
        out_file = abspath(join(orig_dir, base + ext + '.gz'))
    else:
        out_file = abspath(join(getcwd(), base + ext + '.gz'))

    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


def matlab_noddi_processing(caps_directory, num_cores, bStep, name='NoddiMatlab'):
    """
        This is a function to fit NODDI onto the multiple shells diffusion data based on this paper: NODDI: Practical in vivo neurite orientation
        dispersion and density imaging of the human brain, Neuroimage, 2012, Gary Zhang.

        TODO: hard code to find the path to matlab script.
        TODO: for the matlab script, using multiple cores to run one subject has a bug for python wrapper, but ok with one core fitting, improve this in the future
        TODO: deal with how to delete the original nii file in the output folder of NODDI
    Args:
        caps_directory: CAPS directory
        tsv: the tsv file containing the participant_id and session_id
        bStep: the bvalue to round
        working_directory: the path to working_directory
        num_cores: number of cores that fit Noddi model
        name: the name of the pipeline

    Returns:

    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    from nipype.algorithms.misc import Gunzip
    import os
    import clinica.pipelines as clp

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['subject_id_list', 'noddi_preprocessed_dwi', 'noddi_preprocessed_bvec',
                'noddi_preprocessed_bval', 'noddi_preprocessed_mask', 'noddi_toolbox_dir', 'nifti_matlib_dir']),
        name='inputnode')

    capsnode = pe.MapNode(name='capsnode',
                       interface=niu.Function(input_names=['subject_id_list', 'caps_directory'],
                                          output_names=['temp_folder'],
                                          function=make_processing_caps), iterfield=['subject_id_list'])
    capsnode.inputs.caps_directory = caps_directory

    path_to_matscript = os.path.join(os.path.dirname(clp.__path__[0]), 'lib/noddi')


    # output node
    outputnode = pe.Node(niu.IdentityInterface(fields=['fit_icvf', 'fit_isovf', 'fit_od']), name='outputnode')

    # gzip the nii.gz img to nii img
    # Gunzip - unzip functional
    gunzip_dwi = pe.MapNode(Gunzip(), name="gunzipdwi", iterfield=['in_file'])
    gunzip_mask = pe.MapNode(Gunzip(), name="gunzipmask", iterfield=['in_file'])

    # Node to wrap noddi matlab toolbox script.
    nodditoolbox = pe.MapNode(name='nodditoolbox',
                       interface=niu.Function(input_names=['output_dir', 'noddi_img', 'brain_mask', 'roi_mask', 'bval', 'bvec', 'prefix', 'bStep', 'num_cores',
                                                       'path_to_matscript', 'noddi_toolbox_dir', 'nifti_matlib_dir'],
                                          output_names=['fit_icvf', 'fit_isovf', 'fit_od'],
                                          function=runmatlab), iterfield=['output_dir', 'noddi_img', 'brain_mask', 'roi_mask', 'bval', 'bvec', 'prefix'])
    nodditoolbox.inputs.path_to_matscript = path_to_matscript
    nodditoolbox.inputs.num_cores = num_cores
    nodditoolbox.inputs.bStep = bStep

    # zip the result imgs
    zip_icvf = pe.MapNode(name='zip_icvf',
                       interface=niu.Function(input_names=['in_file'],
                                          output_names=['out_file'],
                                          function=compress_nii), iterfield=['in_file'])

    zip_isovf = pe.MapNode(name='zip_isovf',
                       interface=niu.Function(input_names=['in_file'],
                                          output_names=['out_file'],
                                          function=compress_nii), iterfield=['in_file'])

    zip_odi = pe.MapNode(name='zip_odi',
                       interface=niu.Function(input_names=['in_file'],
                                          output_names=['out_file'],
                                          function=compress_nii), iterfield=['in_file'])

    # workflow
    nodditoolbox_wf = pe.Workflow(name=name)
    # unzip
    nodditoolbox_wf.connect(inputnode, 'noddi_preprocessed_dwi', gunzip_dwi, 'in_file')
    nodditoolbox_wf.connect(inputnode, 'noddi_preprocessed_mask', gunzip_mask, 'in_file')
    # fit matlab toolbox
    nodditoolbox_wf.connect(gunzip_dwi, 'out_file', nodditoolbox, 'noddi_img')
    nodditoolbox_wf.connect(gunzip_mask, 'out_file', nodditoolbox, 'brain_mask')
    nodditoolbox_wf.connect(gunzip_mask, 'out_file', nodditoolbox, 'roi_mask')
    nodditoolbox_wf.connect(inputnode, 'noddi_preprocessed_bvec', nodditoolbox, 'bvec')
    nodditoolbox_wf.connect(inputnode, 'noddi_preprocessed_bval', nodditoolbox, 'bval')
    nodditoolbox_wf.connect(inputnode, 'subject_id_list', nodditoolbox, 'prefix')
    # nodditoolbox_wf.connect(inputnode, 'bStep', nodditoolbox, 'bStep')
    nodditoolbox_wf.connect(inputnode, 'noddi_toolbox_dir', nodditoolbox, 'noddi_toolbox_dir')
    nodditoolbox_wf.connect(inputnode, 'nifti_matlib_dir', nodditoolbox, 'nifti_matlib_dir')
    nodditoolbox_wf.connect(inputnode, 'subject_id_list', capsnode, 'subject_id_list')
    nodditoolbox_wf.connect(capsnode, 'temp_folder', nodditoolbox, 'output_dir')

    nodditoolbox_wf.connect(nodditoolbox, 'fit_icvf', zip_icvf, 'in_file')
    nodditoolbox_wf.connect(nodditoolbox, 'fit_isovf', zip_isovf, 'in_file')
    nodditoolbox_wf.connect(nodditoolbox, 'fit_od', zip_odi, 'in_file')
    # output node
    nodditoolbox_wf.connect(zip_icvf, 'out_file', outputnode, 'fit_icvf')
    nodditoolbox_wf.connect(zip_isovf, 'out_file', outputnode, 'fit_isovf')
    nodditoolbox_wf.connect(zip_odi, 'out_file', outputnode, 'fit_od')

    return nodditoolbox_wf


def make_processing_caps(subject_id_list, caps_directory):
    """
    This function is used to make the processing caps
    :param subject_id_list:
    :return:
    """
    import os, errno, tempfile

    processed_path = os.path.join(caps_directory, 'subjects', subject_id_list.split('_')[0], subject_id_list.split('_')[1], 'dwi', 'noddi_based_processing')

    try:
        os.makedirs(processed_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    temp_folder = tempfile.mkdtemp()

    return temp_folder
