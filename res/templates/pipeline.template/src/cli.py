"""

"""

import clinica.engine as ce

class {{ pipeline.class_name }}CLI(ce.CmdParser):

    def define_name(self):
        self._name = '{{ pipeline.command_name }}'

    def define_options(self):
        {% for arg in pipeline.args %}

        {% endfor %}

    def run_pipeline(self, args):
        from {{ pipeline.module_name }} import {{ pipeline.class_name }}

        pipeline = FMRIPreprocessing(self.absolute_path(args.bids_directory),
                                     self.absolute_path(args.caps_directory),
                                     self.absolute_path(args.subjects_sessions_tsv))
        pipeline.parameters = {
            'num_slices'        : args.num_slices,
            'time_repetition'   : args.time_repetition,
            'echo_times'        : args.echo_times,
            'blip_direction'    : args.blip_direction,
            'total_readout_time': args.total_readout_time
        }
        pipeline.base_dir = self.absolute_path(args.working_directory)
        pipeline.write_graph()
        pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_threads})