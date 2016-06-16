from nipype import Workflow

class DecoreWorkflow(Workflow):
    tasks = {}

    def __init__(self, name, base_dir=None):
        Workflow.__init__(self, name, base_dir)

    def add_task(self, name, type):
        if self.tasks.has_key(name):
            self.tasks[name].append(type)
        else:
            self.tasks[name] = [type]

    def run(self, plugin=None, plugin_args=None, updatehash=False):
        Workflow.run(self, plugin, plugin_args, updatehash)
        [task(self) for task in self.tasks.get('run', [])]

def dwork(type, name):
    def decorator(func):
        def wrappper(*args, **kwargs):
            wf = func(*args, **kwargs)
            if not isinstance(wf, Workflow):
                raise Exception('Need Workflow instance')
            wf.__class__ = DecoreWorkflow
            wf.add_task(name, type)
            return wf
        return wrappper
    return decorator

def RunDecorator(type):
    return dwork(type, "run")

def Dump(Workflow):
    import pickle
    with open(Workflow.name + ".pkl", "wb") as dump_file:
        pickle.dump(Workflow, dump_file)
