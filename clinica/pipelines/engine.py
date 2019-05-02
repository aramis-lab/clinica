# coding: utf8

"""

"""

import abc

from nipype.pipeline.engine import Workflow


def get_subject_session_list(input_dir, ss_file=None, is_bids_dir=True, use_session_tsv=False):
    """Parses a BIDS or CAPS directory to get the subjects and sessions.

    This function lists all the subjects and sessions based on the content of
    the BIDS or CAPS directory or (if specified) on the provided
    subject-sessions TSV file.

    Args:
        input_dir: A BIDS or CAPS directory path.
        ss_file: A subjects-sessions file (.tsv format).
        is_bids_dir: Indicates if input_dir is a BIDS or CAPS directory
        use_session_tsv (boolean): Specify if the list uses the sessions listed in the sessions.tsv files

    Returns:
        subjects: A subjects list.
        sessions: A sessions list.
    """
    import clinica.iotools.utils.data_handling as cdh
    import pandas as pd
    import tempfile
    from time import time, strftime, localtime
    import os

    if not ss_file:
        output_dir = tempfile.mkdtemp()
        timestamp = strftime('%Y%m%d_%H%M%S', localtime(time()))
        tsv_file = '%s_subjects_sessions_list.tsv' % timestamp
        ss_file = os.path.join(output_dir, tsv_file)

        cdh.create_subs_sess_list(
            input_dir=input_dir,
            output_dir=output_dir,
            file_name=tsv_file,
            is_bids_dir=is_bids_dir,
            use_session_tsv=use_session_tsv)

    ss_df = pd.io.parsers.read_csv(ss_file, sep='\t')
    if 'participant_id' not in list(ss_df.columns.values):
        raise Exception('No participant_id column in TSV file.')
    if 'session_id' not in list(ss_df.columns.values):
        raise Exception('No session_id column in TSV file.')
    subjects = list(ss_df.participant_id)
    sessions = list(ss_df.session_id)

    # Remove potential whitespace in participant_id or session_id
    return [ses.strip(' ') for ses in sessions], [sub.strip(' ') for sub in subjects]


def postset(attribute, value):
    """Sets the attribute of an object after the execution.

    Args:
        attribute: An object's attribute to be set.
        value: A desired value for the object's attribute.

    Returns:
        A decorator executed after the decorated function is.
    """
    def postset_decorator(func):
        def func_wrapper(self, *args, **kwargs):
            res = func(self, *args, **kwargs)
            setattr(self, attribute, value)
            return res
        return func_wrapper
    return postset_decorator


class Pipeline(Workflow):
    """Clinica Pipeline class.

    This class overwrites the `Workflow` to integrate and encourage the
    use of BIDS and CAPS data structures as inputs and outputs of the pipelines
    developed for the Clinica software.

    The global architecture of a Clinica pipelines is as follow:
        [ Data Input Stream ]
                |
            [ Input ]
                |
            [[ Core ]] <- Could be one or more `npe.Node`s
                |
            [ Output ]
                |
        [ Data Output Stream ]

    Attributes:
        is_built (bool): Informs if the pipelines has been built or not.
        parameters (dict): Parameters of the pipelines.
        info (dict): Information presented in the associated `info.json` file.
        bids_layout (:obj:`BIDSLayout`): Object representing the BIDS directory.
        input_node (:obj:`npe.Node`): Identity interface connecting inputs.
        output_node (:obj:`npe.Node`): Identity interface connecting outputs.
        bids_directory (str): Directory used to read the data from, in BIDS.
        caps_directory (str): Directory used to read/write the data from/to,
            in CAPS.
        subjects (list): List of subjects defined in the `subjects.tsv` file.
            # TODO(@jguillon): Check the subjects-sessions file name.
        sessions (list): List of sessions defined in the `subjects.tsv` file.
        tsv_file (str): Path to the subjects-sessions `.tsv` file.
        info_file (str): Path to the associated `info.json` file.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None):
        """Inits a Pipeline object.

        Args:
            bids_directory (optional): Path to a BIDS directory.
            caps_directory (optional): Path to a CAPS directory.
            tsv_file (optional): Path to a subjects-sessions `.tsv` file.
            name (optional): A pipelines name.
        """
        import inspect
        import os
        self._is_built = False
        self._bids_directory = bids_directory
        self._caps_directory = caps_directory
        self._verbosity = 'debug'
        self._tsv_file = tsv_file
        self._info_file = os.path.join(
            os.path.dirname(os.path.abspath(inspect.getfile(self.__class__))),
            'info.json')
        self._info = {}

        if name:
            self._name = name
        else:
            self._name = self.__class__.__name__
        self._parameters = {}

        if self._bids_directory is None:
            if self._caps_directory is None:
                raise IOError('%s does not contain BIDS nor CAPS directory' %
                              self._name)
            if not os.path.isdir(self._caps_directory):
                raise IOError('The CAPS parameter is not a folder (given path:%s)' %
                              self._caps_directory)

            self._sessions, self._subjects = get_subject_session_list(
                input_dir=self._caps_directory,
                ss_file=self._tsv_file,
                is_bids_dir=False
            )
        else:
            if not os.path.isdir(self._bids_directory):
                raise IOError('The BIDS parameter is not a folder (given path:%s)' %
                              self._bids_directory)
            self._sessions, self._subjects = get_subject_session_list(
                input_dir=self._bids_directory,
                ss_file=self._tsv_file,
                is_bids_dir=True
            )

        self.init_nodes()

    def init_nodes(self):
        """Inits the basic workflow and I/O nodes necessary before build.

        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        if self.get_input_fields():
            self._input_node = npe.Node(name="Input",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=False))
        else:
            self._input_node = None

        if self.get_output_fields():
            self._output_node = npe.Node(name="Output",
                                         interface=nutil.IdentityInterface(
                                             fields=self.get_output_fields(),
                                             mandatory_inputs=False))
        else:
            self._output_node = None

        Workflow.__init__(self, self._name)
        if self.input_node:
            self.add_nodes([self.input_node])
        if self.output_node:
            self.add_nodes([self.output_node])

    def has_input_connections(self):
        """Checks if the Pipeline's input node has been connected.

        Returns:
            True if the input node is connected, False otherwise.
        """
        if self.input_node:
            return self._graph.in_degree(self.input_node) > 0
        else:
            return False

    def has_output_connections(self):
        """Checks if the Pipeline's output node has been connected.

        Returns:
            True if the output node is connected, False otherwise.
        """
        if self.output_node:
            return self._graph.out_degree(self.output_node) > 0
        else:
            return False

    @postset('is_built', True)
    def build(self):
        """Builds the core, input and output nodes of the Pipeline.

        This method first checks it has already been run. It then checks
        the pipelines dependencies and, in this order, builds the core nodes,
        the input node and, finally, the ouput node of the Pipeline.

        Since this method returns the concerned object, it can be chained to
        any other method of the Pipeline class.

        Returns:
            self: A Pipeline object.
        """
        if not self.is_built:
            self.check_dependencies()
            if not self.has_input_connections():
                self.build_input_node()
            self.build_core_nodes()
            if not self.has_output_connections():
                self.build_output_node()
        return self

    def run(self, plugin=None, plugin_args=None, update_hash=False, bypass_check=False):
        """Executes the Pipeline.

        It overwrites the default Workflow method to check if the
        Pipeline is built before running it. If not, it builds it and then
        run it.
        It also checks whether there is enough space left on the disks, and if
        the number of threads to run in parallel is consistent with what is
        possible on the CPU

        Args:
            Similar to those of Workflow.run.

        Returns:
            An execution graph (see Workflow.run).
        """
        if not self.is_built:
            self.build()
        self.check_not_cross_sectional()
        if not bypass_check:
            self.check_size()
            plugin_args = self.update_parallelize_info(plugin_args)
            plugin = 'MultiProc'
        return Workflow.run(self, plugin, plugin_args, update_hash)

    def load_info(self):
        """Loads the associated info.json file.

        Todos:
            - [ ] Raise an appropriate exception when the info file can't open

        Raises:
            None. # TODO(@jguillon)

        Returns:
            self: A Pipeline object.
        """
        import json
        with open(self.info_file) as info_file:
            self.info = json.load(info_file)
        return self

    def check_dependencies(self):
        """Checks if listed dependencies are present.

        Loads the pipelines related `info.json` file and check each one of the
        dependencies listed in the JSON "dependencies" field. Its raises
        exception if a program in the list does not exist or if environment
        variables are not properly defined.

        Todos:
            - [ ] MATLAB toolbox dependency checking
            - [x] check MATLAB
            - [ ] Clinica pipelines dependency checkings
            - [ ] Check dependencies version

        Raises:
            Exception: Raises an exception when bad dependency types given in
            the `info.json` file are detected.

        Returns:
            self: A Pipeline object.
        """
        import clinica.utils.check_dependency as chk

        # Checking functions preparation
        check_software = {
            # 'matlab': chk.check_matlab,
            'ants': chk.check_ants,
            'spm': chk.check_spm,
            'freesurfer': chk.check_freesurfer,
            'fsl': chk.check_fsl,
            'mrtrix': chk.check_mrtrix,
            'matlab': chk.check_matlab
        }
        check_binary = chk.is_binary_present
        # check_toolbox = chk.is_toolbox_present
        # check_pipeline = chk.is_pipeline_present

        # Load the info.json file
        if not self.info:
            self.load_info()

        # Dependencies checking
        for d in self.info['dependencies']:
            if d['type'] == 'software':
                check_software[d['name']]()
            elif d['type'] == 'binary':
                check_binary(d['name'])
            elif d['type'] == 'toolbox':
                pass
            elif d['type'] == 'pipeline':
                pass
            else:
                raise Exception("Unknown dependency type: '%s'." % d['type'])

        self.check_custom_dependencies()

        return self

    def check_size(self):
        """ Checks if the pipeline has enough space on the disk for both
        working directory and caps directory

        Author : Arnaud Marcoux"""
        from os import statvfs
        from os.path import dirname, abspath, join
        from pandas import read_csv
        from clinica.utils.stream import cprint
        from colorama import Fore
        import sys

        SYMBOLS = {
            'customary': ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
            'customary_ext': (
                'byte', 'kilo', 'mega', 'giga', 'tera', 'peta', 'exa',
                'zetta', 'iotta'),
            'iec': ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
            'iec_ext': ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi',
                        'zebi', 'yobi'),
        }

        def bytes2human(n,
                        format='%(value).1f %(symbol)s',
                        symbols='customary'):
            """
            Convert n bytes into a human readable string based on format.
            symbols can be either "customary", "customary_ext", "iec" or "iec_ext",
            see: http://goo.gl/kTQMs

            License :
            Bytes-to-human / human-to-bytes converter.
            Based on: http://goo.gl/kTQMs
            Working with Python 2.x and 3.x.
            Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
            License: MIT
            """
            n = int(n)
            if n < 0:
                raise ValueError("n < 0")
            symbols = SYMBOLS[symbols]
            prefix = {}
            for i, s in enumerate(symbols[1:]):
                prefix[s] = 1 << (i + 1) * 10
            for symbol in reversed(symbols[1:]):
                if n >= prefix[symbol]:
                    value = float(n) / prefix[symbol]
                    return format % locals()
            return format % dict(symbol=symbols[0], value=n)

        def human2bytes(s):
            """
            Attempts to guess the string format based on default symbols
            set and return the corresponding bytes as an integer.
            When unable to recognize the format ValueError is raised.

            License :
            Bytes-to-human / human-to-bytes converter.
            Based on: http://goo.gl/kTQMs
            Working with Python 2.x and 3.x.
            Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
            License: MIT
            """
            init = s
            num = ""
            while s and s[0:1].isdigit() or s[0:1] == '.':
                num += s[0]
                s = s[1:]
            num = float(num)
            letter = s.strip()
            for name, sset in SYMBOLS.items():
                if letter in sset:
                    break
            else:
                if letter == 'k':
                    # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
                    sset = SYMBOLS['customary']
                    letter = letter.upper()
                else:
                    raise ValueError("can't interpret %r" % init)
            prefix = {sset[0]: 1}
            for i, s in enumerate(sset[1:]):
                prefix[s] = 1 << (i + 1) * 10
            return int(num * prefix[letter])

        # Get the number of sessions
        n_sessions = len(self.subjects)
        try:
            caps_stat = statvfs(self.caps_directory)
        except FileNotFoundError:
            # CAPS folder may not exist yet
            caps_stat = statvfs(dirname(self.caps_directory))
        try:
            wd_stat = statvfs(dirname(self.parameters['wd']))
        except (KeyError, FileNotFoundError):
            # Not all pipelines has a 'wd' parameter
            # todo : maybe more relevant to always take base_dir ?
            wd_stat = statvfs(dirname(self.base_dir))

        # Estimate space left on partition/disk/whatever caps and wd is located
        free_space_caps = caps_stat.f_bavail * caps_stat.f_frsize
        free_space_wd = wd_stat.f_bavail * wd_stat.f_frsize

        # space estimation file location
        info_pipelines = read_csv(join(dirname(abspath(__file__)),
                                       'space_required_by_pipeline.csv'),
                                  sep=';')
        pipeline_list = list(info_pipelines.pipeline_name)
        try:
            idx_pipeline = pipeline_list.index(self._name)
            space_needed_caps_1_session = info_pipelines.space_caps[idx_pipeline]
            space_needed_wd_1_session = info_pipelines.space_wd[idx_pipeline]
            space_needed_caps = n_sessions * human2bytes(space_needed_caps_1_session)
            space_needed_wd = n_sessions * human2bytes(space_needed_wd_1_session)
            error = ''
            if free_space_caps == free_space_wd:
                if space_needed_caps + space_needed_wd > free_space_wd:
                    # We assume this is the same disk
                    error = error \
                            + 'Space needed for CAPS and working directory (' \
                            + bytes2human(space_needed_caps + space_needed_wd) \
                            + ') is greater than what is left on your hard drive (' \
                            + bytes2human(free_space_wd) + ')'
            else:
                if space_needed_caps > free_space_caps:
                    error = error + ('Space needed for CAPS (' + bytes2human(space_needed_caps)
                                     + ') is greater than what is left on your hard '
                                     + 'drive (' + bytes2human(free_space_caps) + ')\n')
                if space_needed_wd > free_space_wd:
                    error = error + ('Space needed for working_directory ('
                                     + bytes2human(space_needed_wd) + ') is greater than what is left on your hard '
                                     + 'drive (' + bytes2human(free_space_wd) + ')\n')
            if error != '':
                cprint(Fore.RED + '[SpaceError] ' + error + Fore.RESET)
                while True:
                    cprint('Do you still want to run the pipeline ? (yes/no): ')
                    answer = input('')
                    if answer.lower() in ['yes', 'no']:
                        break
                    else:
                        cprint('Possible answers are yes or no.\n')
                if answer.lower() == 'yes':
                    cprint('Running the pipeline anyway.')
                if answer.lower() == 'no':
                    cprint('Exiting clinica...')
                    sys.exit()

        except ValueError:
            cprint(Fore.RED + 'No info on how much size the pipeline takes. '
                   + 'Running anyway...' + Fore.RESET)

    def update_parallelize_info(self, plugin_args):
        """ Peforms some checks of the number of threads given in parameters,
        given the number of CPUs of the machine in which clinica is running.
        We force the use of plugin MultiProc

        Author : Arnaud Marcoux"""
        from clinica.utils.stream import cprint
        from multiprocessing import cpu_count
        from colorama import Fore

        # count number of CPUs
        n_cpu = cpu_count()
        # Use this var to know in the end if we need to ask the user
        # an other number
        ask_user = False

        try:
            # if no --n_procs arg is used, plugin_arg is None
            # so we need a try / except block
            n_thread_cmdline = plugin_args['n_procs']
            if n_thread_cmdline > n_cpu:
                cprint(Fore.RED + '[Warning] You are trying to run clinica '
                       + 'with a number of threads (' + str(n_thread_cmdline)
                       + ') superior to your number of CPUs (' + str(n_cpu)
                       + ').' + Fore.RESET)
                ask_user = True
        except TypeError:
            cprint(Fore.RED + '[Warning] You did not specify the number of '
                   + 'threads to run in parallel (--n_procs argument).'
                   + Fore.RESET)
            cprint(Fore.RED + 'Computation time can be shorten as you have '
                   + str(n_cpu) + ' CPUs on this computer. We recommand using '
                   + str(n_cpu - 1) + ' threads.\n' + Fore.RESET)
            ask_user = True

        if ask_user:
            while True:
                # While True allows to ask indefinitely until
                # user gives a answer that has the correct format
                # here, positive integer
                cprint('How many threads do you want to use ?')
                answer = input('')
                if answer.isnumeric():
                    if int(answer) > 0:
                        break
                cprint(Fore.RED + 'Your answer must be a positive integer.'
                       + Fore.RESET)
            # If plugin_args is None, create the dictionnary
            # If it already a dict, just update (or create -it is the same
            # code-) the correct key / value
            if plugin_args is None:
                plugin_args = {'n_procs': int(answer)}
            else:
                plugin_args['n_procs'] = int(answer)

        return plugin_args

    def check_not_cross_sectional(self):
        """
        This function checks if the dataset is longitudinal. If it is cross
        sectional, clinica proposes to convert it in a clinica compliant form.

        author: Arnaud Marcoux
        """
        from os import listdir
        from os.path import join, isdir, dirname, abspath, basename
        from colorama import Fore
        from clinica.utils.stream import cprint
        import sys

        def convert_cross_sectional(bids_in,
                                    bids_out,
                                    cross_subjects,
                                    long_subjects):
            """
            This function converts a cross-sectional-bids dataset into a
            longitudinal clinica-compliant dataset

            :param bids_in: cross sectional bids dataset you want to convert
            :param bids_out: converted longitudinal bids dataset
            :param cross_subjects: list of subjects in cross sectional form
            (they need some adjustment)
            :param long_subjects: list of subjects in longitudinal form (they
            just need to be copied)
            :return: nothing
            """
            from os.path import exists, isfile
            from os import mkdir
            from shutil import copytree, copy2

            def add_ses(f):
                """
                Use this function to transform a cross sectional filename into
                a longitudinal one.
                Examples :
                sub-ADNI001_scans.tsv -> sub-ADNI001_ses-M00_scans.tsv
                sub-023a_ses-M12_T1w.nii.gz -> sub-023a_ses-M12_T1w.nii.gz (no
                    modification done if filename already has a session)


                :param f: filename
                :return: filename with '_ses-M00_ added just after
                participant_id
                """
                import re
                # If filename contains ses-..., returns the original filename
                # Regex explication :
                # ^ start of string
                # ([a-zA-Z0-9]*) matches any number of characters from a to z,
                #       A to Z, 0 to 9, and store it in group(1)
                # (?!ses-[a-zA-Z0-9]) do not match if there is already a 'ses-'
                # (.*) catches the rest of the string
                m = re.search(r'(^sub-[a-zA-Z0-9]*)_(?!ses-[a-zA-Z0-9])(.*)',
                              f)
                try:
                    return m.group(1) + '_ses-M00_' + m.group(2)
                except AttributeError:
                    # If something goes wrong, we return the original filename
                    return f

            def copy2_add_ses(src, dst):
                """
                copy2_add_ses calls copy2 function from shutil, but modifies
                the filename of the copied files if they match the regex
                template described in add_ses() function

                :param src: path to the file that needs to be copied
                :param dst: original destination for the copied file
                :return: copy2 with modified filename
                """
                from shutil import copy2
                from os.path import join, dirname, basename

                dst_modified = join(dirname(dst), add_ses(basename(src)))
                return copy2(src, dst_modified)

            if not exists(bids_out):
                # Create the output folder if it does not exists yet
                mkdir(bids_out)

            # First part of the algorithm : deal with subjects that does not
            # have longitudinal (session) information
            for subj in cross_subjects:
                # Get list of of files/folders to copy. Remove hidden element
                # ( if they start with a dot '.')
                to_copy = [f for f in listdir(join(bids_in, subj)) if
                           not f.startswith('.')]
                # Always check that folder are existing, otherwise an exception
                # is raised even though that does not prevent the function from
                # working
                if not exists(join(bids_out, subj)):
                    mkdir(join(bids_out, subj))
                if not exists(join(bids_out, subj, 'ses-M00')):
                    mkdir(join(bids_out, subj, 'ses-M00'))
                for el in to_copy:
                    path_el = join(bids_in, subj, el)
                    if not exists(join(bids_out, subj, 'ses-M00', el)):
                        # If the element to copy is a folder...
                        if isdir(path_el):
                            # Copytree is used with a 'custom' copy function,
                            # to give a correct filemename to the files inside
                            # the copied folders
                            copytree(path_el,
                                     join(bids_out, subj, 'ses-M00',
                                          basename(path_el)),
                                     copy_function=copy2_add_ses)
                        # If the element to copy is a file...
                        elif isfile(path_el):
                            # Modify the filename with add_ses function
                            new_filename_wo_ses = add_ses(el)
                            copy2(path_el,
                                  join(bids_out, subj, 'ses-M00',
                                       new_filename_wo_ses))
            # Second part of the algorithm : deal with subjects that do not
            # have the problem. We only xopy the content of the folder, and no
            # filename needs to be changed
            for su in long_subjects:
                # Do not copy hidden files
                to_copy = [f for f in listdir(join(bids_in, su))
                           if not f.startswith('.')]
                if not exists(join(bids_out, su)):
                    mkdir(join(bids_out, su))
                for el in to_copy:
                    path_el = join(bids_in, su, el)
                    if not exists(join(bids_out, su, el)):
                        # 2 possible cases : element to pcopy is a folder or a
                        # file
                        if isdir(path_el):
                            copytree(path_el,
                                     join(bids_out, su, basename(path_el)))
                        elif isfile(path_el):
                            copy2(path_el,
                                  join(bids_out, su))
        if self.bids_directory is not None:
            bids_dir = abspath(self.bids_directory)
            # Extract all subjects in BIDS directory : element must be a folder and
            # a name starting with 'sub-'
            all_subs = [f for f in listdir(bids_dir)
                        if isdir(join(bids_dir, f)) and f.startswith('sub-')]
            cross_subj = []
            long_subj = []
            for sub in all_subs:
                # Use is_cross_sectional to know if the current subject needs to be
                # labeled as cross sectional
                is_cross_sectional = False
                folder_list = [f for f in listdir(join(bids_dir, sub))
                               if isdir(join(bids_dir, sub, f))]
                for fold in folder_list:
                    if not fold.startswith('ses-'):
                        is_cross_sectional = True
                if is_cross_sectional:
                    cross_subj.append(sub)
                else:
                    long_subj.append(sub)

            # The following code is run if cross sectional subjects have been found
            if len(cross_subj) > 0:
                cprint(Fore.RED + 'It has been determined that '
                       + str(len(cross_subj)) + ' subjects  in your '
                       + 'BIDS folder did not respect the longitudinal '
                       + 'organisation from BIDS specification. Clinica does not '
                       + 'know how to handle cross sectional dataset, but it can '
                       + 'convert it to a Clinica compliant form (using session '
                       + 'ses-M00)\n' + Fore.RESET)
                proposed_bids = join(dirname(bids_dir),
                                     basename(bids_dir) + '_clinica_compliant')

                while True:
                    cprint('Do you want to proceed to the conversion in an other '
                           + 'folder ? (your original BIDS folder will not be'
                           + ' modified, the folder ' + proposed_bids
                           + ' will be created) (yes/no): ')
                    answer = input('')
                    if answer.lower() in ['yes', 'no']:
                        break
                    else:
                        cprint('Possible answers are yes or no.\n')
                if answer.lower() == 'no':
                    cprint('Exiting clinica...')
                    sys.exit()
                else:
                    cprint(
                        'Converting cross-sectional dataset into longitudinal...')
                    convert_cross_sectional(bids_dir,
                                            proposed_bids,
                                            cross_subj,
                                            long_subj)
                    cprint(
                        Fore.GREEN + 'Conversion succeeded. Your clinica-compliant'
                        + ' dataset is located here : ' + proposed_bids
                        + Fore.RESET)

    @property
    def is_built(self): return self._is_built

    @is_built.setter
    def is_built(self, value): self._is_built = value

    @property
    def parameters(self): return self._parameters

    @parameters.setter
    def parameters(self, value):
        self._parameters = value
        # Need to rebuild input, output and core nodes
        self.is_built = False
        self.init_nodes()

    @property
    def info(self): return self._info

    @info.setter
    def info(self, value): self._info = value

    @property
    def bids_layout(self):
        from bids.grabbids import BIDSLayout
        return BIDSLayout(self.bids_directory)

    @property
    def input_node(self): return self._input_node

    @property
    def output_node(self): return self._output_node

    @property
    def bids_directory(self): return self._bids_directory

    @property
    def caps_directory(self): return self._caps_directory

    @property
    def subjects(self): return self._subjects

    @property
    def sessions(self): return self._sessions

    @property
    def tsv_file(self): return self._tsv_file

    @property
    def info_file(self): return self._info_file

    @abc.abstractmethod
    def build_core_nodes(self):
        """Builds the Pipeline's core nodes.

        This method should use and connect to the `Pipeline.input_node` one
        or more core `Node`s. The outputs of the core processing should then be
        connected to the `Pipeline.output_node`.
        """
        pass

    @abc.abstractmethod
    def build_input_node(self):
        """Builds the Pipeline's input data stream node.

        Warnings:
            This method does not modify the `Pipeline.input_node` (see the
            notes about the global architecture in the class documentation).
        """
        pass

    @abc.abstractmethod
    def build_output_node(self):
        """Builds the Pipeline's output data stream node.

        Warnings:
            This method does not modify the `Pipeline.output_node` (see the
            notes about the global architecture in the class documentation).
        """
        pass

    @abc.abstractmethod
    def get_input_fields(self):
        """Lists the input fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.input_node`.
        """
        pass

    @abc.abstractmethod
    def get_output_fields(self):
        """Lists the output fields of the Pipeline.

        Returns:
            A list of strings defining the fields of the `IdentityInterface`
            of the `Pipeline.output_node`.
        """
        pass

    @abc.abstractmethod
    def check_custom_dependencies(self):
        """Checks dependencies provided by the developer.
        """
        pass
