import networkx as nx
import matplotlib.pyplot as plt 
from grave import plot_network
from grave.style import use_attributes
from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import fetch_simulations
import ase, ase.io
from ase.visualize import view
import pymongo, json
from pprint import pprint as pp
import argparse
import numpy as np
#import matplotlib.pyplot as plt 
#import matplotlib
from scipy.spatial.distance import cdist
import sys
from PIL import Image
from collections import OrderedDict as odict

def compute_dg_diff(idx, branch_dct, ne_dct, 
    reference_energy = 0.0, free_energy_correction = 0.0):
    """Helper function to compute the adsorption free energy

    In the case of hydrogen the property is:
    dGdiff(nH) = G(nH) - G(nH-1)) - 0.5 G(H2)

    Args:
        idx (int) : simulation id
        branch_dct (dict) : parent simulation : list of child simulations
        ne_dct (dict) : stores total energies of all calculations with respect to
                        the number of adsorbates and their ids
        reference_energy (float) :  reference energy for the adsorbate. Can be the
                                    total energy of the isolated adsorbate molecule
                                    or a different reference point
        free_energy_correction (float) :    free energy correction of the adsorption 
                                            reaction at hand

        Returns:
        tuple : adsorption free energy (float), 
                assignment of property to 'parent' or 'child' simulation
    """
    # search where idx appears in branch_dct
    for key, ids in zip(branch_dct.keys(), branch_dct.values()):
        if int(idx) in list(ids):
            parent_idx = key
            break
    # get energies     
    for n_adsorbates, v in ne_dct.items():
        if str(idx) in v.keys():
            energy = v[str(idx)]
            child_n_adsorbates = int(n_adsorbates)
        if str(parent_idx) in v.keys():
            parent_energy = v[str(parent_idx)]
            parent_n_adsorbates = int(n_adsorbates)
            min_energy = np.min(list(v.values()))
            parent_energy = min_energy

    assert (parent_n_adsorbates - child_n_adsorbates) == 1 or (parent_n_adsorbates - child_n_adsorbates) == -1

    if parent_n_adsorbates > child_n_adsorbates:
        # dg diff
        dg_diff = parent_energy - energy - reference_energy + free_energy_correction 
        assignment = "parent"
    else:
        dg_diff = energy - parent_energy - reference_energy + free_energy_correction 
        assignment = "child"
    print("#################################################")
    print("dG_diff", dg_diff, assignment)
    print("dG_diff components", energy, parent_energy, reference_energy ,free_energy_correction )
    print("n_adsorbates", "parent", parent_n_adsorbates, "child", child_n_adsorbates)
    print("#################################################")
    return dg_diff, assignment


def get_external_database(**extdb_connect):
    #extdb_connect["username"] = username
    #extdb_connect["password"] = password
    extdb_connect["host"] = extdb_connect.get("host",
        "nanolayers.dyndns.org:27017")
    extdb_connect["db_name"] = extdb_connect.get("db_name",
        "testdb")        
    extdb_connect["authsource"] = extdb_connect.get("authsource",
        extdb_connect["db_name"])
    CLIENT = pymongo.MongoClient(extdb_connect["host"], 
            username = extdb_connect["username"],
            password = extdb_connect["password"],
            authSource = extdb_connect["authsource"])
    db = CLIENT[extdb_connect["db_name"]]
    return db, extdb_connect

def plot_atoms(atoms, rotation, filename, is_render_skipped = False):
    from ase import Atoms
    from ase.io import write
    from ase.data.colors import jmol_colors

    cell = atoms.get_cell()
    cell[:] = 0
    atoms.set_cell(cell)

    #Make colors
    #colors = np.zeros((len(atoms), 3))
    #elements = atoms.get_atomic_numbers()[cluster_ids]
    #colors[cluster_ids, :] = jmol_colors[elements, :]

    # Textures
    tex = ['jmol',] * 288 + ['glass',] * 288+ ['ase3',] * 288 + ['vmd',] * 288 

    # keywords
    kwargs = { # Keywords that exist for eps, png, and pov
    'rotation': rotation,
    'show_unit_cell': 0,
    'colors': None,
    'radii': None,
    }

    extra_kwargs = { # For povray files only
    'display'      : False, # Display while rendering
    'pause'        : False, # Pause when done rendering (only if display)
    'transparent'  : True, # Transparent background
    'canvas_width' : 500,  # Width of canvas in pixels
    'canvas_height': None,  # Height of canvas in pixels
    'camera_dist'  : 50.,   # Distance from camera to front atom
    'image_plane'  : None,  # Distance from front atom to image plane
                            # (focal depth for perspective)
    'camera_type'  : 'perspective', # perspective, ultra_wide_angle
    'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
    'area_light'   : [(2., 3., 40.) ,# location
                      'White',       # color
                      .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
    #'background'   : 'White',        # color
    'textures'     : tex, # Length of atoms list of texture names
    'celllinewidth': 0.05, # Radius of the cylinders representing the cell
    }
    # Make flat png file
    write(filename + '_flat.png', atoms, **kwargs)
    
    kwargs.update(extra_kwargs)
    # Make the raytraced image
    if is_render_skipped == True:
        #print("skipping image rendering")
        return filename + '_flat.png'
    else:
        #print("rendering image:", filename)
        write(filename + '.pov', atoms, run_povray=True, **kwargs)
        return filename + '.png'


def load_image(filename):
    import numpy as np
    import matplotlib.pyplot as plt 
    from scipy.misc import imread, imsave
    image_data = imread(filename).astype(np.float32)
    print('Image Size: ', image_data.size)
    print('Image Shape: ', image_data.shape)

    scaled_image_data = image_data / 255.
    #imsave('test_out.png', scaled_image_data)
    #plt.imshow(scaled_image_data)
    #plt.show()
    return scaled_image_data


def concatenate_images_horizontally(filenames):
    images = [Image.open(x) for x in filenames]
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
      new_im.paste(im, (x_offset,0))
      x_offset += im.size[0]

    new_im.save('test.jpg')
    return new_im

def get_external_database(**extdb_connect):
    #extdb_connect["username"] = username
    #extdb_connect["password"] = password
    extdb_connect["host"] = extdb_connect.get("host",
        "nanolayers.dyndns.org:27017")
    extdb_connect["db_name"] = extdb_connect.get("db_name",
        "testdb")        
    extdb_connect["authsource"] = extdb_connect.get("authsource",
        extdb_connect["db_name"])
    CLIENT = pymongo.MongoClient(extdb_connect["host"], 
            username = extdb_connect["username"],
            password = extdb_connect["password"],
            authSource = extdb_connect["authsource"])
    db = CLIENT[extdb_connect["db_name"]]
    return db, extdb_connect


def hilighter(event):
    # if we did not hit a node, bail
    if not hasattr(event, 'nodes') or not event.nodes:
        return

    # pull out the graph,
    graph = event.artist.graph

    # clear any non-default color on nodes
    for node, attributes in graph.nodes.data():
        attributes.pop('color', None)

    for u, v, attributes in graph.edges.data():
        attributes.pop('width', None)

    for node in event.nodes:
        graph.nodes[node]['color'] = 'C1'
        #print(node)

        for edge_attribute in graph[node].values():
            edge_attribute['width'] = 3 

    # update the screen
    event.artist.stale = True
    event.artist.figure.canvas.draw_idle()
    print(event.nodes)

    simulation_id = event.nodes[0]
    simulation = simulations[str(simulation_id)]
    pp(simulation)
    print("simulation id", simulation_id)
    atoms = atoms_dict_to_ase(simulation["atoms"])
    view(atoms)
    

def draw_workflow(dol, simulations):
    graph = nx.from_dict_of_lists(dol, create_using=None)
    
    fig, ax = plt.subplots()
    art = plot_network(graph, ax=ax, node_style=use_attributes(),
            edge_style=use_attributes(), 
            #layout="kamada_kawai",
            #node_label_style = {'font_weight': 'bold', 'font_size': 15}
            )
    art.set_picker(10)
    ax.set_title('workflow simulations, click on the nodes to see structure')
    fig.canvas.mpl_connect('pick_event', hilighter)
    plt.show()
    plt.close()
    return 


def get_simulations_from_workflow(workflow_id, is_test):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db, _ = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    cursor = ext_db["simulations"].find({"workflow_id": workflow_id},)
    simulations_list = [document for document in cursor]
    simulation_ids = [str(document["_id"]) for document in simulations_list]
    simulations  = dict(zip(simulation_ids, simulations_list))
    return simulations


def get_dol_from_simulations(simulations, is_test = False):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"
    ext_db, _ = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    dol = {}
    new_simulations = simulations.copy()
    for key, simulation in simulations.items():
        source_id = simulation["source_id"]
        dol[key] = [str(source_id)]
        try:
            source_simulation = simulations[str(source_id)]
        except:
            print("getting source simulation from database")
            source_simulation  = ext_db["simulations"].find_one({"_id": source_id})
            new_simulations[str(source_id)] = source_simulation

    return dol, new_simulations


def reverse(dct):
    """
    Return a reversed dict, with common values in the original dict
    grouped into a list in the returned dict.

    Example:
    >>> d = ReversibleDict({'a': 3, 'c': 2, 'b': 2, 'e': 3, 'd': 1, 'f': 2})
    >>> d.reversed()
    {1: ['d'], 2: ['c', 'b', 'f'], 3: ['a', 'e']}
    """

    revdict = {}
    for k, v in dct.items():
        for entry in v:
            revdict.setdefault(entry, []).append(k)
    return revdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--draw', action='store_true')
    parser.add_argument('--skiprender', action='store_true')
    parser.add_argument('-j','--jsons', nargs='+', help='list of Fireworks json files', default = [])

    args = parser.parse_args()
    print(args.ids)
    is_test = args.test
    is_draw = args.draw
    is_render_skipped = args.skiprender

    workflow_id = args.ids[0]

    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"
        db, _ = get_external_database(
            username = "myusername",
            password = "mypassword",
            db_name = db_name,
            )

    """wf_output = {"branch_dct" : branch_dct,
         "property" : energies,
         "calc_ids" : calc_ids,
         "ne_dct" : ne_dct,
         "calc_parents" : calc_parents,
         "root_history" : fw_spec["temp"]["root_history"],
         "step_history" : fw_spec["temp"]["step_history"],
         "free_energy_correction" : fw_spec["temp"]["free_energy_correction"],
         "reference_energy" : fw_spec["temp"]["reference_energy"],
         "direction" : direction,
         "is_return" : is_return,
         "n_adsorbates" : n_adsorbates,
        }
    """

    workflows = db['workflows']
    #returned = workflows.update_one({'_id':workflow_id}, {"$set": {"output" : wf_output}}, upsert=False)
    dct = workflows.find_one({'_id': workflow_id})
    #pp(dct)
    dct = dct["output"]
    step_history = dct["step_history"]
    root_history = dct["root_history"]
    branch_dct = dct["branch_dct"]
    ne_dct = dct["ne_dct"]
    reference_energy = dct["reference_energy"]
    free_energy_correction = dct["free_energy_correction"]

    # prepare energy plot
    x_n_adsorbates = []
    y_dg = []
    x_epochs = []
    x_root = []

    # prepare plotting grid

    MAX_VERTICAL = 60
    MAX_HORIZONTAL = 7
    fig, ax = plt.subplots(MAX_VERTICAL,MAX_HORIZONTAL, gridspec_kw=dict(wspace=0.0, hspace=0.1), figsize=(MAX_HORIZONTAL * 2, MAX_VERTICAL * 1.1))
    for i in range(MAX_VERTICAL):
        for k in range(MAX_HORIZONTAL):
            ax[i, k].set_xticks([])
            ax[i, k].set_yticks([])

    simulations = {}    
    for idx in args.ids:
        new_simulations = get_simulations_from_workflow(idx, is_test)
        simulations.update(new_simulations)
    
    dol, simulations = get_dol_from_simulations(simulations)
    print("dictionary of links")
    print(dol)
    if is_draw:
        draw_workflow(dol, simulations)

    # find root
    reverse_dol = reverse(dol)
    print("REVERSE dictionary of links")
    print(reverse_dol)
    root_dol = reverse_dol.copy()

    children_lst = []
    for key, lst in reverse_dol.items():
        children_lst.extend(lst)

    for key, lst in reverse_dol.items():
        if key in children_lst:
            root_dol.pop(key)

    print(root_dol.keys())
    if len(root_dol.keys()) > 1:
        print("WARNING. More than one root simulation. Proceeding with first entry")

    root_id = list(root_dol.keys())[0]
    print(root_id)

    children_ids = [root_id]
    is_root = True

    altered_atoms = {}

    i = 0
    epoch = 0
    while i < MAX_VERTICAL:
        current_simulations = odict()
        print("# simulations on current level", len(children_ids))
        print(children_ids)
        for idx in children_ids:
            current_simulations[idx] = simulations[idx]

        current_structures = odict()
        for idx, simulation in current_simulations.items():
            step = simulation["operations"][0]
            #print("step:", step)

            atoms = atoms_dict_to_ase(simulation["atoms"])
            current_structures[idx] = atoms
            n_adsorbates = len(atoms[atoms.get_atomic_numbers() != 0]) - 55
            #print(n_adsorbates)


        # compare parent with current structures
        marked_structures = odict()
        if not is_root: # add safeguard for root and reprints

            for idx, atoms in current_structures.items():
                # mapping of parents to children
                parent_idx = current_simulations[idx]["source_id"]
                step = current_simulations[idx]["operations"][0]

                parent_atoms = parent_structures[str(parent_idx)]

                dmat = cdist(atoms.get_positions(), parent_atoms.get_positions())
                #print(dmat.min(axis = 0))
                if step == {'add_adsorbate': -1}:

                    removed_atom_id = np.argmax(dmat.min(axis = 0))
                    marked_structure = parent_atoms.copy()
                    marked_structure.symbols[removed_atom_id] = 'X'
                    marked_structures[idx] = marked_structure

                elif step == {'add_adsorbate': 1}:
                    added_atom_id = np.argmax(dmat.min(axis = 1))
                    marked_structure = atoms.copy()
                    marked_structure.symbols[added_atom_id] = 'Be'
                    marked_structures[idx] = marked_structure
                    altered_atoms[str(idx)] = added_atom_id

                    #view(marked_structure)
                    #view(parent_atoms)
                    #view(atoms)
                    #exit(1)
                elif (step == 'cp2k'):
                    #print("ENTERING step cp2k")
                    added_atom_id = altered_atoms.get(str(parent_idx), None)
                    marked_structure = atoms.copy()
                    if added_atom_id:
                        marked_structure.symbols[int(added_atom_id)] = 'Be'
                    marked_structures[idx] = marked_structure

                    print(idx)

                    try:
                        dg, assignment = compute_dg_diff(idx, branch_dct, ne_dct, 
                            reference_energy = reference_energy, free_energy_correction = free_energy_correction)
                        print(dg, assignment)
                        n_adsorbates = len(marked_structure) - 55
                        if assignment == "parent":
                            # dGdiff = dG(n) - dG(n-1) - u
                            n_adsorbates += 1
                        if (dg > -100) and (dg < 100):
                            x_epochs.append(epoch)
                            x_n_adsorbates.append(n_adsorbates)
                            y_dg.append(dg * 27.2114)
                            if int(idx) in root_history:
                                x_root.append(1)
                            else:
                                x_root.append(0)

                    except:
                        print("skipping ", idx)
                else:
                    marked_structures =  current_structures
        else:
            marked_structures =  current_structures

        for j, (simulation_id, atoms) in enumerate(marked_structures.items()):
            filename_root = "ladder" + "_coverage"
            filename1 = plot_atoms(atoms, rotation = '', 
                filename = filename_root, is_render_skipped = is_render_skipped)
            filename2 = plot_atoms(atoms, rotation = '180y', 
                filename = filename_root + "_rot", is_render_skipped = is_render_skipped)
            current_image = concatenate_images_horizontally([filename1, filename2])
            n_adsorbates = len(atoms[atoms.get_atomic_numbers() != 0]) - 55

            k = j % MAX_HORIZONTAL
            if (j % MAX_HORIZONTAL == 0) and (j != 0):
                i += 1

            if i >= MAX_VERTICAL:
                print("WARNING: cutting current series of images")
                pass
            else:   
                im = ax[i, k].imshow(current_image, interpolation='none',
                    origin='lower',
                    #extent=[-2, 4, -3, 2], clip_on=True
                    )
                #ax[i, k].set_title(str(i) + " " + str(j))
                if int(simulation_id) in root_history:
                    color = 'r'
                else:
                    color = 'b'
                ax[i,k].text(i, k, str(n_adsorbates) + "   " + str(simulation_id), fontsize=12,
                       ha="left", va="center", color=color)
                #ax[i,k].annotate('figure fraction', xy=(.1, .5), xycoords='figure fraction',
                    #horizontalalignment='left', verticalalignment='top',
                    #fontsize=12)


        # save parent structures

        parent_simulations = current_simulations
        parent_structures = current_structures

        # get children

        children_ids = []
        #print(reverse_dol)
        #print(current_simulations.keys())
        for idx, _ in current_simulations.items():
            children_ids.extend(reverse_dol.get(str(idx), []))
            #print("parent, children", idx, reverse_dol.get(str(idx), []))

        if len(list(set(children_ids))) != len(children_ids):
            print("removing duplicate children")
            children_ids = list(set(children_ids))
        #print(children_ids)

        is_root = False
        i += 1
        epoch += 1
        print("EPOCH", epoch)
    
    # finalize
    plt.savefig("overview_coverage_ladder.png", dpi = 300)
    #plt.show()

    plt.close()

    x_root = np.array(x_root, dtype=int)
    #print("root datapoints:", x_root.sum())
    #print("root datapoints:", x_root)#

    #print(np.array(x_epochs)[x_root])
    #print(np.array(y_dg)[x_root])

    #print(np.array(x_epochs))
    #print(np.array(y_dg))

    plt.scatter(x_n_adsorbates, y_dg)
    plt.scatter(np.array(x_n_adsorbates)[x_root == 1], np.array(y_dg)[x_root == 1])
    plt.xlabel('# adsorbates')
    plt.ylabel(r'$\Delta \ G \ [eV]$')
    plt.tight_layout()
    plt.savefig("energies_coverage_ladder.png", dpi = 300)
    plt.show()
    plt.close()

    #import matplotlib.cm as cm
    #colors = cm.rainbow(x_n_adsorbates)
    cm = plt.cm.get_cmap('RdYlBu')

    sc = plt.scatter(x_epochs, y_dg, c = x_n_adsorbates, cmap = cm)
    plt.scatter(np.array(x_epochs)[x_root == 1], np.array(y_dg)[x_root == 1],s=80, facecolors='r', edgecolors='r')
    plt.xlabel('epoch count')
    plt.ylabel(r'$\Delta \ G \ [eV]$')
    plt.colorbar(sc)
    plt.tight_layout()
    plt.savefig("history_energies_coverage_ladder.png", dpi = 300)
    plt.show()
