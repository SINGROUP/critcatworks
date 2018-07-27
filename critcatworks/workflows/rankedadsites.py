from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib

# internal modules
from critcatworks.clusgeo import get_adsites
from critcatworks.database import read_structures

@explicit_serialize
class SimpleTestTask(FiretaskBase):
   """ 
   Simple FireTask to see what attributes FireTasks have.

   Args:
       base_name (str): A random base name.
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

@explicit_serialize
class DFTSetupTask(FiretaskBase):
   """ 
   Task to setup DFT calculations.

   Args:
        None
   """

   _fw_name = 'DFTSetupTask'
   required_params = []
   optional_params = []

   def run_task(self, fw_spec):
       #print(self.__dict__)
       print("DFTSetupTask not implemented yet")
       pp(fw_spec)
       print("forwarding information")
       pass_spec = fw_spec
       print("dummy outputs generated")

       print("TODO getting structures and energies")

       mod_spec =[{'_set' : {'dft_outputs->id' + str(dft_params[0]) : dft_params[0]}}]
       return FWAction(update_spec =pass_spec, mod_spec=mod_spec)


def branch_dft_calculations_fw():
    task = DFTSetupTask(spec={'dft_params':[]}, inputs=['dft_params'])
    task_dict = task.to_dict()
    firetask1 = ForeachTask(task=task_dict, split='dft_params', spec={'dft_params':[1,2,3]})
    fw = Firework([firetask1])
    return fw


def test_foreachtask_workflow():
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
    # create the Firework consisting of multiple tasks
    firetask1 = TemplateWriterTask({'context': {'opt1': 5.0, 'opt2': 'fast method'}, 'template_file': 'simple_template.txt', 'output_file': 'inputs.txt'})
    firetask2 = ScriptTask.from_str('wc -w < inputs.txt > words.txt')
    firetask3 = FileTransferTask({'files': [{'src': 'words.txt', 'dest': '~/words.txt'}], 'mode': 'copy'})
    wf = Firework([firetask1, firetask2, firetask3])
    return wf

def cp2k_test_workflow():
    simple_task = SimpleTestTask(spec={'simple': 'stuff'})
    firetask = PyTask(run_cp2k("cp2k_mm_energy.inp"))
    wf = Firework([simple_task, firetask])

    return wf


def run_cp2k(template_file):
    #export ASE_CP2K_COMMAND="mpirun -n 2 cp2k_shell.popt"
    from ase.calculators.cp2k import CP2K
    from ase.build import molecule

    with open('cp2k_mm_energy.inp', 'r') as f:
        content = f.read()
    calc = CP2K(inp=content)
    atoms = molecule('H2O', calculator=calc)
    atoms.center(vacuum=2.0)
    print(atoms.get_potential_energy())
    return

def get_adsites_workflow(path, adsorbate_energy=0.0, adsorbate_name='H'):
    """
    Workflow to determine the adsorption sites and energies of a set of
    nanocluster structures using CP2K and Clusgeo
    """
    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information
    abspath = pathlib.Path(path).resolve()
    fw_read_structures = read_structures(abspath)
    # FireWork: Determine adsites and add to database
    fw_get_adsites = get_adsites(adsorbate_energy=0.0, adsorbate_name='H')
    # FireWork: FPS ranking

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)

    # FireWork: update database, 
    # (includes reading relaxed structure and energy)

    # FireWork: machine learning from database

    # FireWork: check if converged, give intermediary overview.
    # give summary when finished


    wf = Workflow([fw_read_structures, fw_get_adsites], links_dict = {fw_read_structures: [fw_get_adsites]})
    return wf




if __name__ == "__main__":
    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)
    #launchpad = LaunchPad(host="myhost", port=12345, \
    #name="fireworks_testing_db", username="my_user", \
    #password="my_pass")

    #wf = get_adsites_workflow()
    #wf = dummy_workflow()
    #wf = cp2k_test_workflow()
    #wf = test_foreachtask_workflow()
    wf = get_adsites_workflow(path = "/l/programs/critcatworks/tests/dummy_db/nc_structures/")

    # store workflow and launch it locally, single shot
    launchpad.add_wf(wf)


    # excecute workflow
    IS_QUEUE = False
    if IS_QUEUE:
        launch_rocket_to_queue(launchpad, FWorker(), adapter, launcher_dir=mypath, reserve=True)
    else:
        #launch_rocket(launchpad, FWorker())
        rapidfire(launchpad, FWorker())


    if IS_QUEUE:
        # recover offline fireworks
        for i in range(0,10):
            time.sleep(5)
            ids =launchpad.get_fw_ids()
            for idx in ids:
                launchpad.recover_offline(launch_id = idx)



