"""
Script that generate the skeleton for a new Pipeline
"""
from clinica.engine.cmdparser import CmdParser
from os.path import join, isdir, splitext
from os import mkdir, getcwd


class CmdGenerateTemplates(CmdParser):
    """
    Using the jinja2 library and the pipelines template's file (in the resource dir)
    the user can create his own run-to-go pipelines
    """

    def __init__(self):
        CmdParser.__init__(self)

    def define_name(self):
        self._name = 'template'

    def define_description(self):
        self._description = 'Generate the skeleton for a new pipeline (for developers)'

    def define_options(self):
        self._args.add_argument("name",
                                help='The pipelines title')
        self._args.add_argument("-d", "--output_dir",
                                help='Define the path where generate the directory')

    def run_command(self, args):
        # Parsing input arguments
        if args.output_dir is None:
            args.output_dir = getcwd()
        pipeline = dict()
        pipeline['title'] = args.name
        pipeline['module_name'] = pipeline['title'].replace(' ', '_').lower()
        pipeline['class_name'] = pipeline['title'].replace(' ', '')
        pipeline['command_name'] = pipeline['title'].replace(' ', '-').lower()
        pipeline['dir'] = join(args.output_dir, pipeline['module_name'])

        from jinja2 import Environment, PackageLoader, select_autoescape
        env = Environment(
            loader=PackageLoader('clinica', 'resources/templates/pipeline_template'),
            autoescape=select_autoescape(['py'])
        )

        if not isdir(pipeline['dir']):
            mkdir(pipeline['dir'])

        for template_file in env.list_templates(extensions='j2'):
            template = env.get_template(template_file)
            rendered_file, template_ext = splitext(template_file)
            rendered_filename, rendered_file_ext = splitext(rendered_file)
            if rendered_file_ext == '.py' and rendered_filename != '__init__':
                path_to_write = join(pipeline['dir'], pipeline['module_name'] + '_' + rendered_file)
            else:
                path_to_write = join(pipeline['dir'], rendered_file)
            with open(path_to_write, 'w+') as file_to_write:
                print("Generating template %s" % file_to_write.name)
                file_to_write.write(template.render(pipeline=pipeline))
