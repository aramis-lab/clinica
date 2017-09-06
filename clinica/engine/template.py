"""
Prototype in order to generate automatic Pipeline
"""
from clinica.engine.cmdparser import CmdParser
import abc
import string
from os.path import realpath,split,join,isdir,splitext
from os import mkdir, getcwd, listdir


def to_camel_case(text):
    return ''.join(x for x in text.title() if not x.isspace())


class AbstractTemplate:
    __metaclass__ = abc.ABCMeta
    args = None

    def path_template(self):
        return join(split(realpath(__file__))[0], '../resources/template')

    @abc.abstractmethod
    def file_template(self): pass

    @abc.abstractmethod
    def file_to_generate(self, args): pass

    def add_argument(self, args): pass

    @abc.abstractmethod
    def dictionary(self, args): pass


class PipelineTemplate(AbstractTemplate):

    def file_template(self):
        return 'pipeline.template'

    def file_to_generate(self, args):
        return '%sPipeline.py' % args.name

    def dictionary(self, args):
        return {
            'clinica_major_version': 0,
            'clinica_minor_version': 1,
            'clinica_patch_version': 0,
            'pipeline_name': args.name
        }


class VisualizeTemplate(AbstractTemplate):

    def file_template(self):
        return 'visualize.template'

    def file_to_generate(self, args):
        return '%sVisualize.py' % args.name

    def dictionary(self, args):
        return {
            'clinica_major_version': 0,
            'clinica_minor_version': 1,
            'clinica_patch_version': 0,
            'pipeline_name': args.name
        }


class CmdGenerateTemplates(CmdParser):

    def __init__(self):
        CmdParser.__init__(self)

    def define_name(self):
        self._name = 'template'

    def define_options(self):
        self._args.add_argument("name",
                                help='The pipeline title')
        self._args.add_argument("-d", "--output_dir",
                                help='Define the path where generate the directory')

    def run_pipeline(self, args):
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
                print "Generating template %s" % file_to_write.name
                file_to_write.write(template.render(pipeline=pipeline))
