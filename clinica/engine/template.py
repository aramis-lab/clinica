"""
Proptotype in order to generate automatique Pipeline
"""
from clinica.engine.cmdparser import CmdParser
import abc
import string
from os.path import realpath,split,join,isdir
from os import mkdir, getcwd


def CamelCase(text):
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
        return {'clinica_major_version' : 1,
                'clinica_minor_version' : 0,
                'clinica_patch_version' : 0,
                'pipeline_name': args.name
                }

class VisualizeTemplate(AbstractTemplate):

    def file_template(self):
        return 'visualize.template'

    def file_to_generate(self, args):
        return '%sVisualize.py' % args.name

    def dictionary(self, args):
        return {'clinica_major_version' : 1,
                'clinica_minor_version' : 0,
                'clinica_patch_version' : 0,
                'pipeline_name': args.name
                }

class CmdGenerateTemplates(CmdParser):

    def __init__(self):
        CmdParser.__init__(self)

    def define_name(self):
        self._name = 'template'

    def define_options(self):
        self.template = [PipelineTemplate(), VisualizeTemplate()]
        for template in self.template : template.add_argument(self._args)
        self._args.add_argument("--name",
                               required=True,
                               help='The pipeline name')
        self._args.add_argument("--generate_dir",
                                help='Define the path where generate the directory')

    def run_pipeline(self, args):
        args.name = args.name.strip()

        if args.generate_dir is None:
            args.generate_dir = join(getcwd(), 'Clinica%s' % args.name)

        if isdir(args.generate_dir) is False:
            mkdir(args.generate_dir)

        for template in self.template :
            with open(join(template.path_template(),template.file_template()),'r') as file_to_read:
                with open(join(args.generate_dir, template.file_to_generate(args)), 'w') as file_to_write:
                    print "Generating template %s" % file_to_write.name
                    file_to_write.write(string.Template(file_to_read.read()).substitute(template.dictionary(args)))


