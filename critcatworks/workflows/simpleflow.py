# this workflow is for testing if fireworks runs
from fireworks import Firework, FWorker, LaunchPad, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask

import os,time
from pprint import pprint as pp
import pathlib
from critcatworks.database import mylaunchpad

@explicit_serialize
class SimpleTestTask(FiretaskBase):
   """ 
   Simple FireTask to see what attributes FireTasks have.
   """

   _fw_name = 'SimpleTestTask'
   required_params = []
   optional_params = []

   def run_task(self, fw_spec):
       #print(self.__dict__)
       print("This is a simple test task")
       #print("outputs", self.outputs)
       pp(fw_spec)
       #ranking = fw_spec["ranking"]
       #ranking.append(1)
       dft_params = fw_spec["dft_params"]
       #dft_params.append(1)
       update_spec ={'dft_params': dft_params}
       mod_spec =[{'_set' : {'dft_outputs->id' + str(dft_params[0]) : dft_params[0]}}]
       #return FWAction(update_spec=update_spec)
       return FWAction(update_spec ={'dft_params': dft_params}, mod_spec=mod_spec)


def test_foreachtask_workflow():
    """
    Workflow to test fireworks ForeachTask with a dummy workflow
    """
    simple_task = SimpleTestTask(spec={'dft_params':[]}, inputs=['dft_params'])
    simple_dict = simple_task.to_dict()
    #resimple_task = SimpleTestTask.from_dict(simple_dict)
    #print(simple_task)
    print(simple_dict)
    #print(resimple_task)
    firetask1 = ForeachTask(task=simple_dict, split='dft_params', spec={'dft_params':[1,2,3]})

    fw1 = Firework([simple_task], spec={'dft_params':[-1,-2,-3], 'dft_outputs':{}   }, fw_id=1)
    fw2 = Firework([firetask1],   fw_id=2)
    fw3 = Firework([firetask1],   fw_id=3)
    fw4 = Firework([simple_task], fw_id=4)

    #wf = Firework([simple_task, firetask1, simple_task], spec={'dft_params':[1,2,3], 'ranking': [-3,-1] })
    workflow = Workflow([fw1, fw2, fw4], {1: [2], 2: [4],})
    return workflow

def dummy_workflow():
    """
    dummy fireworks Workflow
    """
    # create the Firework consisting of multiple tasks
    firetask1 = TemplateWriterTask({'context': {'opt1': 5.0, 'opt2': 'fast method'}, 'template_file': 'simple_template.txt', 'output_file': 'inputs.txt'})
    firetask2 = ScriptTask.from_str('wc -w < inputs.txt > words.txt')
    firetask3 = FileTransferTask({'files': [{'src': 'words.txt', 'dest': '~/words.txt'}], 'mode': 'copy'})
    wf = Firework([firetask1, firetask2, firetask3])
    return wf

if __name__ == "__main__":
    import logging
    IS_QUEUE = False
    print("dummy workflow")
    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logdir = str(pathlib.Path(".").resolve())
        logging.basicConfig(filename = logdir + "/logfile_dummy_workflow.log", level=logging.INFO)

    # set up the LaunchPad and reset it
    launchpad = mylaunchpad.create_launchpad()
    launchpad.reset('', require_password=False)
    wf = dummy_workflow()

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)


    #launch_rocket(launchpad, FWorker())
    rapidfire(launchpad, FWorker())
