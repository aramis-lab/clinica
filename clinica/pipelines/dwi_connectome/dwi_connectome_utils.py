# coding: utf8

# The 3 imports are required for @jguillon's merge requests.

from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base
from nipype.interfaces.base import traits, TraitedSpec, File, Undefined

import os.path as op


def get_luts():
    import os
    from clinica.utils.exceptions import ClinicaException

    try:
        # For aparc+aseg.mgz file:
        default = os.path.join(os.environ['FREESURFER_HOME'],
                               'FreeSurferColorLUT.txt')
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(os.environ['FREESURFER_HOME'],
                              'FreeSurferColorLUT.txt')

        # TODO: Add custom Lausanne2008 LUTs here.
    except KeyError:
        raise ClinicaException('Could not find FREESURFER_HOME environment variable.')

    return [default, a2009s]


def get_conversion_luts():
    import os
    from clinica.utils.exceptions import ClinicaException

    try:
        # For aparc+aseg.mgz file:
        default = os.path.join(os.environ['MRTRIX_HOME'],
                               'share/mrtrix3/labelconvert/fs_default.txt')
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(os.environ['MRTRIX_HOME'],
                              'share/mrtrix3/labelconvert/fs_a2009s.txt')

        # TODO: Add custom Lausanne2008 conversion LUTs here.
    except KeyError:
        raise ClinicaException('Could not find MRTRIX_HOME environment variable.')

    return [default, a2009s]


def get_containers(subjects, sessions):

    return [
        'subjects/' + subjects[i] + '/' + sessions[i] + '/dwi'
        for i in range(len(subjects))
    ]


def get_caps_filenames(dwi_file):

    import re

    m = re.search(r'\/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_preproc', dwi_file)

    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')

    source_file = m.group(1)

    response = source_file + '_model-CSD_responseFunction.txt'
    fod = source_file + '_model-CSD_diffmodel.mif'
    tracts = source_file + '_model-CSD_tractography.tck'
    nodes = [source_file + '_atlas-desikan_parcellation.nii.gz',
             source_file + '_atlas-destrieux_parcellation.nii.gz']
    # TODO: Add custom Lausanne2008 node files here.
    connectomes = [source_file + '_model-CSD_atlas-desikan_connectivity.tsv',
                   source_file + '_model-CSD_atlas-destrieux_connectivity.tsv']
    # TODO: Add custom Lausanne2008 connectome files here.

    return response, fod, tracts, nodes, connectomes


def print_begin_pipeline(in_bids_or_caps_file):
    """
    """
    from clinica.utils.stream import cprint
    import re
    import datetime
    from colorama import Fore

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)',
                  in_bids_or_caps_file)
    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.')
    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s Running pipeline for %s...' % (
        Fore.BLUE, now, Fore.RESET, m.group(0)))


def print_end_pipeline(in_bids_or_caps_file, final_file):
    """
    """
    from clinica.utils.stream import cprint
    import re
    import datetime
    from colorama import Fore

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)',
                  in_bids_or_caps_file)
    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.')
    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s ...%s has completed.' % (
        Fore.GREEN, now, Fore.RESET, m.group(0)))


# From this line, everything that is defined was to avoid waiting for nipype's
# team to merge @jguillon's requests in the `master` branch and also to avoid
# having to update clinica with the last version of nipype.

class EstimateFODInputSpec(MRTrix3BaseInputSpec):
    algorithm = traits.Enum(
        'csd',
        'msmt_csd',
        argstr='%s',
        position=-8,
        mandatory=True,
        desc='FOD algorithm')
    in_file = File(
        exists=True,
        argstr='%s',
        position=-7,
        mandatory=True,
        desc='input DWI image')
    wm_txt = File(
        argstr='%s', position=-6, mandatory=True, desc='WM response text file')
    wm_odf = File(
        'wm.mif',
        argstr='%s',
        position=-5,
        usedefault=True,
        mandatory=True,
        desc='output WM ODF')
    gm_txt = File(argstr='%s', position=-4, desc='GM response text file')
    gm_odf = File('gm.mif', usedefault=False, argstr='%s',
                  position=-3, desc='output GM ODF')
    csf_txt = File(argstr='%s', position=-2, desc='CSF response text file')
    csf_odf = File('csf.mif', usedefault=False, argstr='%s',
                   position=-1, desc='output CSF ODF')
    mask_file = File(exists=True, argstr='-mask %s', desc='mask image')

    # DW Shell selection options
    shell = traits.List(
        traits.Float,
        sep=',',
        argstr='-shell %s',
        desc='specify one or more dw gradient shells')
    max_sh = traits.Int(
        8, usedefault=True,
        argstr='-lmax %d',
        desc='maximum harmonic degree of response function')
    in_dirs = File(
        exists=True,
        argstr='-directions %s',
        desc=('specify the directions over which to apply the non-negativity '
              'constraint (by default, the built-in 300 direction set is '
              'used). These should be supplied as a text file containing the '
              '[ az el ] pairs for the directions.'))


class EstimateFODOutputSpec(TraitedSpec):
    wm_odf = File(argstr='%s', desc='output WM ODF')
    gm_odf = File(argstr='%s', desc='output GM ODF')
    csf_odf = File(argstr='%s', desc='output CSF ODF')


class EstimateFOD(MRTrix3Base):
    """
    Estimate fibre orientation distributions from diffusion data using spherical deconvolution

    Example
    -------

    >>> import nipype.interfaces.mrtrix3 as mrt
    >>> fod = mrt.EstimateFOD()
    >>> fod.inputs.algorithm = 'csd'
    >>> fod.inputs.in_file = 'dwi.mif'
    >>> fod.inputs.wm_txt = 'wm.txt'
    >>> fod.inputs.grad_fsl = ('bvecs', 'bvals')
    >>> fod.cmdline                               # doctest: +ELLIPSIS
    'dwi2fod -fslgrad bvecs bvals -lmax 8 csd dwi.mif wm.txt wm.mif gm.mif csf.mif'
    >>> fod.run()                                 # doctest: +SKIP
    """

    _cmd = 'dwi2fod'
    input_spec = EstimateFODInputSpec
    output_spec = EstimateFODOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['wm_odf'] = op.abspath(self.inputs.wm_odf)
        if self.inputs.gm_odf != Undefined:
            outputs['gm_odf'] = op.abspath(self.inputs.gm_odf)
        if self.inputs.csf_odf != Undefined:
            outputs['csf_odf'] = op.abspath(self.inputs.csf_odf)
        return outputs


class TractographyInputSpec(MRTrix3BaseInputSpec):
    sph_trait = traits.Tuple(
        traits.Float,
        traits.Float,
        traits.Float,
        traits.Float,
        argstr='%f,%f,%f,%f')

    in_file = File(
        exists=True,
        argstr='%s',
        mandatory=True,
        position=-2,
        desc='input file to be processed')

    out_file = File(
        'tracked.tck',
        argstr='%s',
        mandatory=True,
        position=-1,
        usedefault=True,
        desc='output file containing tracks')

    algorithm = traits.Enum(
        'iFOD2',
        'FACT',
        'iFOD1',
        'Nulldist',
        'SD_Stream',
        'Tensor_Det',
        'Tensor_Prob',
        usedefault=True,
        argstr='-algorithm %s',
        desc='tractography algorithm to be used')

    # ROIs processing options
    roi_incl = traits.Either(
        File(exists=True),
        sph_trait,
        argstr='-include %s',
        desc=('specify an inclusion region of interest, streamlines must'
              ' traverse ALL inclusion regions to be accepted'))
    roi_excl = traits.Either(
        File(exists=True),
        sph_trait,
        argstr='-exclude %s',
        desc=('specify an exclusion region of interest, streamlines that'
              ' enter ANY exclude region will be discarded'))
    roi_mask = traits.Either(
        File(exists=True),
        sph_trait,
        argstr='-mask %s',
        desc=('specify a masking region of interest. If defined,'
              'streamlines exiting the mask will be truncated'))

    # Streamlines tractography options
    step_size = traits.Float(
        argstr='-step %f',
        desc=('set the step size of the algorithm in mm (default is 0.1'
              ' x voxelsize; for iFOD2: 0.5 x voxelsize)'))
    angle = traits.Float(
        argstr='-angle %f',
        desc=('set the maximum angle between successive steps (default '
              'is 90deg x stepsize / voxelsize)'))
    n_tracks = traits.Int(
        argstr='-select %d',
        max_ver=0.4,
        desc=('set the desired number of tracks. The program will continue'
              ' to generate tracks until this number of tracks have been '
              'selected and written to the output file'))
    select = traits.Int(
        argstr='-select %d',
        min_ver=3,
        desc=('set the desired number of tracks. The program will continue'
              ' to generate tracks until this number of tracks have been '
              'selected and written to the output file'))
    max_tracks = traits.Int(
        argstr='-maxnum %d',
        desc=('set the maximum number of tracks to generate. The program '
              'will not generate more tracks than this number, even if '
              'the desired number of tracks hasn\'t yet been reached '
              '(default is 100 x number)'))
    max_length = traits.Float(
        argstr='-maxlength %f',
        desc=('set the maximum length of any track in mm (default is '
              '100 x voxelsize)'))
    min_length = traits.Float(
        argstr='-minlength %f',
        desc=('set the minimum length of any track in mm (default is '
              '5 x voxelsize)'))
    cutoff = traits.Float(
        argstr='-cutoff %f',
        desc=('set the FA or FOD amplitude cutoff for terminating '
              'tracks (default is 0.1)'))
    cutoff_init = traits.Float(
        argstr='-initcutoff %f',
        desc=('set the minimum FA or FOD amplitude for initiating '
              'tracks (default is the same as the normal cutoff)'))
    n_trials = traits.Int(
        argstr='-trials %d',
        desc=('set the maximum number of sampling trials at each point'
              ' (only used for probabilistic tracking)'))
    unidirectional = traits.Bool(
        argstr='-unidirectional',
        desc=('track from the seed point in one direction only '
              '(default is to track in both directions)'))
    init_dir = traits.Tuple(
        traits.Float,
        traits.Float,
        traits.Float,
        argstr='-initdirection %f,%f,%f',
        desc=('specify an initial direction for the tracking (this '
              'should be supplied as a vector of 3 comma-separated values'))
    noprecompt = traits.Bool(
        argstr='-noprecomputed',
        desc=('do NOT pre-compute legendre polynomial values. Warning: this '
              'will slow down the algorithm by a factor of approximately 4'))
    power = traits.Int(
        argstr='-power %d',
        desc=('raise the FOD to the power specified (default is 1/nsamples)'))
    n_samples = traits.Int(
        4, usedefault=True,
        argstr='-samples %d',
        desc=('set the number of FOD samples to take per step for the 2nd '
              'order (iFOD2) method'))
    use_rk4 = traits.Bool(
        argstr='-rk4',
        desc=('use 4th-order Runge-Kutta integration (slower, but eliminates'
              ' curvature overshoot in 1st-order deterministic methods)'))
    stop = traits.Bool(
        argstr='-stop',
        desc=('stop propagating a streamline once it has traversed all '
              'include regions'))
    downsample = traits.Float(
        argstr='-downsample %f',
        desc='downsample the generated streamlines to reduce output file size')

    # Anatomically-Constrained Tractography options
    act_file = File(
        exists=True,
        argstr='-act %s',
        desc=('use the Anatomically-Constrained Tractography framework during'
              ' tracking; provided image must be in the 5TT '
              '(five - tissue - type) format'))
    backtrack = traits.Bool(
        argstr='-backtrack', desc='allow tracks to be truncated')

    crop_at_gmwmi = traits.Bool(
        argstr='-crop_at_gmwmi',
        desc=('crop streamline endpoints more '
              'precisely as they cross the GM-WM interface'))

    # Tractography seeding options
    seed_sphere = traits.Tuple(
        traits.Float,
        traits.Float,
        traits.Float,
        traits.Float,
        argstr='-seed_sphere %f,%f,%f,%f',
        desc='spherical seed')
    seed_image = File(
        exists=True,
        argstr='-seed_image %s',
        desc='seed streamlines entirely at random within mask')
    seed_rnd_voxel = traits.Tuple(
        File(exists=True),
        traits.Int(),
        argstr='-seed_random_per_voxel %s %d',
        xor=['seed_image', 'seed_grid_voxel'],
        desc=('seed a fixed number of streamlines per voxel in a mask '
              'image; random placement of seeds in each voxel'))
    seed_grid_voxel = traits.Tuple(
        File(exists=True),
        traits.Int(),
        argstr='-seed_grid_per_voxel %s %d',
        xor=['seed_image', 'seed_rnd_voxel'],
        desc=('seed a fixed number of streamlines per voxel in a mask '
              'image; place seeds on a 3D mesh grid (grid_size argument '
              'is per axis; so a grid_size of 3 results in 27 seeds per'
              ' voxel)'))
    seed_rejection = File(
        exists=True,
        argstr='-seed_rejection %s',
        desc=('seed from an image using rejection sampling (higher '
              'values = more probable to seed from'))
    seed_gmwmi = File(
        exists=True,
        argstr='-seed_gmwmi %s',
        requires=['act_file'],
        desc=('seed from the grey matter - white matter interface (only '
              'valid if using ACT framework)'))
    seed_dynamic = File(
        exists=True,
        argstr='-seed_dynamic %s',
        desc=('determine seed points dynamically using the SIFT model '
              '(must not provide any other seeding mechanism). Note that'
              ' while this seeding mechanism improves the distribution of'
              ' reconstructed streamlines density, it should NOT be used '
              'as a substitute for the SIFT method itself.'))
    max_seed_attempts = traits.Int(
        argstr='-max_seed_attempts %d',
        desc=('set the maximum number of times that the tracking '
              'algorithm should attempt to find an appropriate tracking'
              ' direction from a given seed point'))
    out_seeds = File(
        'out_seeds.nii.gz', usedefault=True,
        argstr='-output_seeds %s',
        desc=('output the seed location of all successful streamlines to'
              ' a file'))


class TractographyOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='the output filtered tracks')
    out_seeds = File(
        desc=('output the seed location of all successful'
              ' streamlines to a file'))


class Tractography(MRTrix3Base):
    """
    Performs streamlines tractography after selecting the appropriate
    algorithm.

    .. [FACT] Mori, S.; Crain, B. J.; Chacko, V. P. & van Zijl,
      P. C. M. Three-dimensional tracking of axonal projections in the
      brain by magnetic resonance imaging. Annals of Neurology, 1999,
      45, 265-269

    .. [iFOD1] Tournier, J.-D.; Calamante, F. & Connelly, A. MRtrix:
      Diffusion tractography in crossing fiber regions. Int. J. Imaging
      Syst. Technol., 2012, 22, 53-66

    .. [iFOD2] Tournier, J.-D.; Calamante, F. & Connelly, A. Improved
      probabilistic streamlines tractography by 2nd order integration
      over fibre orientation distributions. Proceedings of the
      International Society for Magnetic Resonance in Medicine, 2010, 1670

    .. [Nulldist] Morris, D. M.; Embleton, K. V. & Parker, G. J.
      Probabilistic fibre tracking: Differentiation of connections from
      chance events. NeuroImage, 2008, 42, 1329-1339

    .. [Tensor_Det] Basser, P. J.; Pajevic, S.; Pierpaoli, C.; Duda, J.
      and Aldroubi, A. In vivo fiber tractography using DT-MRI data.
      Magnetic Resonance in Medicine, 2000, 44, 625-632

    .. [Tensor_Prob] Jones, D. Tractography Gone Wild: Probabilistic Fibre
      Tracking Using the Wild Bootstrap With Diffusion Tensor MRI. IEEE
      Transactions on Medical Imaging, 2008, 27, 1268-1274


    Example
    -------

    >>> import nipype.interfaces.mrtrix3 as mrt
    >>> tk = mrt.Tractography()
    >>> tk.inputs.in_file = 'fods.mif'
    >>> tk.inputs.roi_mask = 'mask.nii.gz'
    >>> tk.inputs.seed_sphere = (80, 100, 70, 10)
    >>> tk.cmdline                               # doctest: +ELLIPSIS
    'tckgen -algorithm iFOD2 -samples 4 -output_seeds out_seeds.nii.gz \
-mask mask.nii.gz -seed_sphere \
80.000000,100.000000,70.000000,10.000000 fods.mif tracked.tck'
    >>> tk.run()                                 # doctest: +SKIP
    """

    _cmd = 'tckgen'
    input_spec = TractographyInputSpec
    output_spec = TractographyOutputSpec

    def _format_arg(self, name, trait_spec, value):
        if 'roi_' in name and isinstance(value, tuple):
            value = ['%f' % v for v in value]
            return trait_spec.argstr % ','.join(value)

        return super(Tractography, self)._format_arg(name, trait_spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs
