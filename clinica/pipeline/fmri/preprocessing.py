"""

"""


import nipype.pipeline.engine as npe
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.utility as nutil
import clinica.pipeline.engine as cpe
from clinica.utils.io import zip_nii, unzip_nii


def filter_data(file_list):
    new_file_list = filter(None, file_list)
    return new_file_list


class FMRIPreprocessing(cpe.Pipeline):
    """Create fMRI preprocessing pipeline object.

    TODO: Long description.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A nipype workflow object containing the full fMRI preprocessing pipeline.

    Raises:
        IOError:
    """

    @cpe.postset('is_built', True)
    def build(self):

        # Data grabbing
        # ================
        dg_node = npe.Node(name='DataGrabbing',
                           interface=nio.DataGrabber(infields=['sessions', 'subjects'],
                                                     outfields=['img', 'hdr']))
        dg_node.inputs.base_directory=self.input_dir
        dg_node.inputs.template='*'
        dg_node.inputs.raise_on_empty=False
        dg_node.inputs.sort_filelist=True
        dg_node.inputs.field_template=dict(img='%s/%s/func/%s_%s_bold.nii.gz',
                                           hdr='%s/%s/func/%s_%s_bold.json')
        dg_node.inputs.template_args=dict(img=[['subjects', 'sessions', 'subjects', 'sessions']],
                                          hdr=[['subjects', 'sessions', 'subjects', 'sessions']])
        dg_node.inputs.sessions = ['ses-M00', 'ses-M00'] # self.sessions
        dg_node.inputs.subjects = ['sub-HMTC2010_07_15_MEMEP_Sujet20',
                                   'sub-HMTC2011_05_13_MEMEP_PAT29'] #self.subjects
        # NOTE: This is going to work only when nipype team will have solved this issue
        #   (https://github.com/nipy/nipype/issues/1783), I removed l.1183-1184.

        # Data filtering
        # =============
        filter_node = npe.Node(name='DataFiltering',
                               interface=nutil.Function(input_names=['file_list'],
                                                        output_names=['new_file_list'],
                                                        function=filter_data))
        # Unzipping
        # =========
        unzip_node = npe.MapNode(name='Unzipping',
                                 iterfield=['in_file'],
                                 interface=nutil.Function(input_names=['in_file'],
                                                          output_names=['out_file'],
                                                          function=unzip_nii))

        # Slice timing correction
        # =======================
        st_node = npe.Node(name="SliceTimingCorrection",
                           interface=spm.SliceTiming())
        st_node.inputs.num_slices = self.parameters['num_slices']
        st_node.inputs.time_repetition = self.parameters['time_repetition']
        st_node.inputs.time_acquisition = self.parameters['time_repetition'] - self.parameters['time_repetition'] \
                                                                               / float(self.parameters['num_slices'])
        st_node.inputs.slice_order = range(1, self.parameters['num_slices']+1)
        st_node.inputs.ref_slice = self.parameters['num_slices'] / 2

        # Motion correction
        # =================
        mc_node = npe.Node(name="MotionCorrection",
                           interface=spm.Realign())
        mc_node.inputs.register_to_mean = True
        mc_node.inputs.write_mask = True

        # Registration
        # ============

        # Normalization
        # =============

        # Smoothing
        # =========

        # Zipping
        # =======
        zip_node = npe.MapNode(name='Zipping',
                               iterfield=['in_file'],
                               interface=nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=zip_nii))

        # Connection
        # ==========
        self.connect([
            (dg_node,       filter_node,    [('img',                    'file_list')]),
            (filter_node,   unzip_node,     [('new_file_list',          'in_file')]),
            # (dg_node,       unzip_node, [('img', 'in_file')]),
            (unzip_node,    st_node,        [('out_file',               'in_files')]),
            (st_node,       mc_node,        [('timecorrected_files',    'in_files')]),
            (mc_node,       zip_node,       [('realigned_files',        'in_file')]),
        ])

        return self
