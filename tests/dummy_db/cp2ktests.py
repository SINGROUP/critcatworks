from fireworks import Firework, FWorker, LaunchPad, ScriptTask, TemplateWriterTask, FileTransferTask
from fireworks.core.rocket_launcher import launch_rocket
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
import os,time

def parse_input(inp):
    """Parses the given CP2K input string"""
    root_section = 'CP2K_INPUT'
    section_stack = [root_section]

    for line in inp.split('\n'):
        line = line.split('!', 1)[0].strip()
        if len(line) == 0:
            continue

        if line.upper().startswith('&END'):
            s = section_stack.pop()
        elif line[0] == '&':
            parts = line.split(' ', 1)
            name = parts[0][1:]
            if len(parts) > 1:
                s = InputSection(name=name, params=parts[1].strip())
            else:
                s = InputSection(name=name)
            section_stack[-1].subsections.append(s)
            section_stack.append(s)
        else:
            section_stack[-1].keywords.append(line)

    return root_section


def old_parse_input(inp):
    """Parses the given CP2K input string"""
    root_section = InputSection('CP2K_INPUT')
    section_stack = [root_section]

    for line in inp.split('\n'):
        line = line.split('!', 1)[0].strip()
        if len(line) == 0:
            continue

        if line.upper().startswith('&END'):
            s = section_stack.pop()
        elif line[0] == '&':
            parts = line.split(' ', 1)
            name = parts[0][1:]
            if len(parts) > 1:
                s = InputSection(name=name, params=parts[1].strip())
            else:
                s = InputSection(name=name)
            section_stack[-1].subsections.append(s)
            section_stack.append(s)
        else:
            section_stack[-1].keywords.append(line)

    return root_section

def dummy_workflow():
    # create the Firework consisting of multiple tasks
    firetask1 = TemplateWriterTask({'context': {'opt1': 5.0, 'opt2': 'fast method'}, 'template_file': 'simple_template.txt', 'output_file': 'inputs.txt'})
    firetask2 = ScriptTask.from_str('wc -w < inputs.txt > words.txt')
    firetask3 = FileTransferTask({'files': [{'src': 'words.txt', 'dest': '~/words.txt'}], 'mode': 'copy'})
    fw = Firework([firetask1, firetask2, firetask3])
    return fw

def run_cp2k(template_file):
    from ase.calculators.cp2k import CP2K, parse_input
    from ase.build import molecule
    from ase.optimize import BFGS

    with open('cp2k_mm_energy.inp', 'r') as f:
        content = f.read()
    calc = CP2K()
    #print("input parameters")
    #print(calc.inp)
    print(calc)
    print("input parameters")
    #print(calc.todict())

    print("parsing input")
    parsed = parse_input(content)
    print(parsed)
    print(parsed.name)
    print(parsed.params)
    print(parsed.keywords)
    print(parsed.subsections)

    exit(0)
    atoms = molecule('H2O', calculator=calc)
    atoms.center(vacuum=2.0)
    gopt = BFGS(atoms, logfile=None)
    gopt.run(fmax=1e-6)

    print(atoms.get_potential_energy())
    return
 

def get_adsites_workflow():
    """
    Workflow to determine the adsorption sites and energies of a set of
    nanocluster structures using CP2K and Clusgeo
    """
    # FireWork: Read nanocluster structures and initialise a database
    # object containing set information

    # FireWork: Determine adsites and add to database

    # FireWork: FPS ranking

    # FireWork: setup, run and extract DFT calculation
    # (involves checking for errors in DFT and rerunning)

    # FireWork: update database, 
    # (includes reading relaxed structure and energy)

    # FireWork: machine learning from database

    # FireWork: check if converged, give intermediary overview.
    # give summary when finished






if __name__ == "__main__":
    run_cp2k('cp2k_mm_energy.inp')

