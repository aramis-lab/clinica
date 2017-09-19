"""
Utility class to serialize data and run pipeline
"""
from nipype import Workflow

class ClinicaWorkflow(Workflow):
    """
    Execute pipeline and store data in order to run the pipeline
    """
    tasks = {}
    is_reset = False

    def __init__(self, name, base_dir=None):
        Workflow.__init__(self, name, base_dir)
        reset(self)

    def reset(self):
        self.data = {}
        self.timestamp = None
        self.dump = False

    def add_task(self, name, type):
        if self.tasks.has_key(name):
            self.tasks[name].append(type)
        else:
            self.tasks[name] = [type]

    def add_data(self, name, info):
        if self.data.has_key(name):
            self.data[name].append(info)
        else:
            self.data[name] = info

    def run(self, plugin=None, plugin_args=None, updatehash=False):
        Workflow.run(self, plugin, plugin_args, updatehash)
        [task(self) for task in self.tasks.get('run', [])]


def dwork(type, name):
    """
    Decorator : program specific task for a pipeline
    """
    def decorator(func):
        def wrappper(*args, **kwargs):
            wf = func(*args, **kwargs)
            if not isinstance(wf, Workflow):
                raise Exception('Need Workflow instance')
            wf.__class__ = ClinicaWorkflow
            if wf.is_reset == False:
                wf.reset()
                wf.is_reset = True
            wf.add_task(name, type)
            return wf
        return wrappper
    return decorator


def RunDecorator(type):
    """
    The RUN decorator
    """
    return dwork(type, "run")


def Dump(clinicaWorkflow):
    """
    The DUMP decorator
    use:
        @Dump
        def my_pipeline(...):
    """
    import cPickle
    import datetime
    from os.path import realpath,join
    print(clinicaWorkflow.base_dir)
    with open(join(realpath(clinicaWorkflow.base_dir), "clinica.pkl"), "wb") as dump_file:
        clinicaWorkflow.timestamp = datetime.datetime.utcnow()
        cPickle.dump(clinicaWorkflow, dump_file)
    return clinicaWorkflow


def Visualize(application, parameters, matches):
    """
    The Visualize decorator
    use:
        @Visualize("freeview", "-t %1")
        def my_pipeline(...):
    """
    def visu(clinicaWorkflow):
        clinicaWorkflow.add_data('visualize', [application, parameters, matches])
        return Dump(clinicaWorkflow)
    return RunDecorator(visu)
