from fireworks import Firework, FWorker, LaunchPad, PyTask, ScriptTask, TemplateWriterTask, FileTransferTask, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.firetasks.dataflow_tasks import ForeachTask
from pprint import pprint as pp
import pathlib

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
