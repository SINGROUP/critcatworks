# utility script for singlesites or molsinglesites
# input three arguments: workflow_id, element1, element2
# Given two element types (in atomic numbers), adsorption energy
# distribution is plotted w.r.t composition
# interactive, inspect structures upon clicking

from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import fetch_simulations
import ase, ase.io
from ase.visualize import view
import pymongo
from pprint import pprint as pp
import argparse
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib import pyplot as plt
import matplotlib.colors
from scipy.stats import gaussian_kde
from ase.data import chemical_symbols as pte
import matplotlib.gridspec as gridspec

def plot_atoms_properties(atoms, rotation, filename, cluster_ids = np.arange(55), minimum = 0.0, maximum = 1.0, is_render_skipped = False):
    from ase import Atoms
    from ase.io import write
    from ase.data.colors import jmol_colors

    cell = atoms.get_cell()
    cell[:] = 0
    atoms.set_cell(cell)

    #Make colors
    colors = np.zeros((len(atoms), 3))
    elements = atoms.get_atomic_numbers()[cluster_ids]
    colors[cluster_ids, :] = jmol_colors[elements, :]
    non_cluster_ids = np.setdiff1d(np.arange(len(atoms)), cluster_ids)

    cmap = plt.cm.rainbow
    norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

    colors[non_cluster_ids] = cmap(norm(atoms.get_masses()[non_cluster_ids]))[:,:3]

    # Textures
    tex = ['jmol',] * 288 + ['glass',] * 288+ ['ase3',] * 288 + ['vmd',] * 288 

    # keywords
    kwargs = { # Keywords that exist for eps, png, and pov
    'rotation': rotation,
    'show_unit_cell': 0,
    'colors': colors,
    'radii': None,
    }

    extra_kwargs = { # For povray files only
    'display'      : False, # Display while rendering
    'pause'        : False, # Pause when done rendering (only if display)
    'transparent'  : True, # Transparent background
    'canvas_width' : 1500,  # Width of canvas in pixels
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
        print("skipping image rendering")
        return filename + '_flat.png'
    else:
        print("rendering image:", filename)
        write(filename + '.pov', atoms, run_povray=True, **kwargs)
        # make colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        fig, ax = plt.subplots()
        sm.set_array([])
        fig.colorbar(sm)
        ax.remove()
        plt.savefig('global_colorbar.png',bbox_inches='tight')
        #plt.show()
        plt.close()
        return filename + '.png'

def interactive_plot_scatter_mean(x, y, x_mean, y_mean, label = "", name = "noname_scatter.png"):
    print("interactive plotting of composition vs adsorption energy")
    # picker
    plt.close('all')
    fig, (ax0) = plt.subplots(1, 1)
    fig.canvas.mpl_connect('pick_event', onpick3)

    xy = np.vstack([x,y])
    #z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    #idx = z.argsort()
    #x, y, z = x[idx], y[idx], z[idx]

    plt.ylabel(r'$\Delta \ E \ [eV]$')
    plt.xlabel(label)
    plt.scatter(x, y, alpha = 0.5, picker = True)

    #plt.plot(x_mean, y_mean, label = "mean " + label, 
    #    color='r', marker = '_', markersize = 20.0)
    
    # bounding box not working
    #plt.legend(bbox_to_anchor=(1.20,1), loc="upper left")
    #plt.legend(loc='upper left')

    #plt.tight_layout()
    #plt.savefig(name, bbox_inches="tight")
    #plt.savefig(name)
    plt.show()
    #plt.clf()
    return

def plot_scatter_mean(x, y, x_mean, y_mean, label = "", name = "noname_scatter.png"):
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    plt.scatter(x, y, c=z, alpha = 0.5)
    plt.ylabel(r'$\Delta \ E \ [eV]$')
    plt.xlabel(label)
    plt.plot(x_mean, y_mean, 
        color='r', marker = '_', markersize = 20.0)
    
    # bounding box not working
    #plt.legend(bbox_to_anchor=(1.20,1), loc="upper left")
    #plt.legend(loc='upper left')

    plt.tight_layout()
    #plt.savefig(name, bbox_inches="tight")
    plt.savefig(name)
    #plt.show()
    plt.clf()
    return None


def get_external_database(**extdb_connect):
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
    return db

def get_simulations_from_elements(el1, el2, workflow_id = 45, is_test = False):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    # subset tests for metal elements plus common adsorbate light elements
    cursor = ext_db.simulations.find(
        {'workflow_id' : workflow_id, 
        "$expr":{"$setIsSubset":["$atoms.numbers",[el1,el2, 0,1,6,7,8,9]]}
        }, 
        {"atoms.numbers" : 1, "output.total_energy" : 1, 
            "output.is_converged" : 1, "source_id" : 1, "operations" : 1,
            "nanoclusters" : 1
            }
        )
    simulations_list = [document for document in cursor]
    simulation_ids = [str(document["_id"]) for document in simulations_list]
    simulations  = dict(zip(simulation_ids, simulations_list))
    print("get simulations from elements")
    return simulations, simulations_list


def get_ml(workflow_id, is_test):
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    ext_db = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    ml = ext_db["machine_learning"].find({"workflow_id": workflow_id}, {"metrics_validation" : 1,"metrics_training" : 1,"metrics_test" : 1,"workflow_id" : 1, 
        "training_set": 1, "test_set" : 1, "prediction_set" : 1, "output" : 1,
        })
    return ml

def fetch_simulations(simulation_ids, is_test = False):
    """Fetches simulation records by simulation id

    Args:
        simulation_ids (1D ndarray) : unique identifiers of the simulation collection.
        is_test (bool)  : switch to use either test or production database

    Returns:
        list : documents of the simulation collection
    """
    if is_test:
        db_name = "testdb"
    else:
        db_name = "ncdb"

    db = get_external_database(
        username = "myusername",
        password = "mypassword",
        db_name = db_name,
        )
    cursor =   db['simulations'].find({"_id" : {"$in" : simulation_ids }}) 
    simulations_list = [document for document in cursor]
    temp_simulation_ids = [str(document["_id"]) for document in simulations_list]
    simulations  = dict(zip(temp_simulation_ids, simulations_list))
    simulations_list = [simulations[str(idx)] for idx in simulation_ids]
    return simulations, simulations_list


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

def onpick3(event):
    ind = event.ind
    print('onpick3 scatter:', ind)
    try:
        first_id = ind[0]
    except:
        return
    print("taking first point", first_id)
    simulation = initial_guess_simulations_list[first_id]
    pp(simulation)
    print("simulation id", simulation["_id"])
    atoms = atoms_dict_to_ase(simulation["atoms"])
    view(atoms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id, element1, element2')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--skiprender', action='store_true')
    parser.add_argument('--noninter', action='store_true')


    # input:    wf id
    #           el1, el2
    args = parser.parse_args()
    print(args.ids)
    is_test = args.test
    is_render_skipped = args.skiprender
    is_noninter = args.noninter
    workflow_id = args.ids[0]
    element1 = int(args.ids[1])
    element2 = int(args.ids[2])
    print("workflow_id", workflow_id)
    print("element 1", element1, pte[element1])
    print("element 2", element2, pte[element2])


    # get labels from last prediction
    ml = get_ml(workflow_id, is_test)
    for entry in ml:
        train_size = len(entry["training_set"])
        test_size = len(entry["test_set"])
        print("train size", train_size, "test_size", test_size)
        #pp(entry)

    print("getting property from last machine learning run")
    test_ids = entry["test_set"]
    training_ids = entry["training_set"]
    prediction_ids = entry["prediction_set"]
    all_ids = []
    all_ids.extend(training_ids)
    all_ids.extend(test_ids)
    all_ids.extend(prediction_ids)

    test_property = entry["output"]["label_test"]
    training_property = entry["output"]["label_train"]
    predicted_property = entry["output"]["label_predicted"]
    y_property = []
    y_property.extend(training_property)
    y_property.extend(test_property)
    y_property.extend(predicted_property)
    y_property = np.array(y_property) * 27.2114 # in eV

    # make dictionary of id : label
    property_dct = {all_ids[i] : y for i, y in enumerate(y_property)}

    # get cp2k simulations containing two types of elements + pure elements #
    #########################################################################
    simulations, simulations_list = get_simulations_from_elements(element1, element2, 
        workflow_id = workflow_id, is_test = is_test)

    # removing not required simulations #
    #####################################
    new_simulations_list = []
    for simulation in simulations_list:
        if simulation["_id"] not in all_ids:
            #print(simulation["_id"], "not in required ids")
            continue
        else:
            new_simulations_list.append(simulation)

    simulations_list = new_simulations_list

    # get reference ids of nanoclusters
    reference_id_dct = {}
    for simulation in simulations_list:
        # get reference id of nanoclusters
        ref_id = simulation["nanoclusters"][0]["reference_id"]
        value = reference_id_dct.get(ref_id, [])
        value.append(int(simulation["_id"]))
        reference_id_dct[ref_id] = value

    print("downloaded", len(simulations), "simulations")


    # get source simulations if necessary
    target_operations = [{'add_adsorbate' : 1}]
    initial_guess_ids = []

    for simulation in simulations_list:
        operations = simulation["operations"]
        if operations == target_operations:
            initial_guess_ids.append(int(simulation["_id"]))
        else:
            source_id = simulation["source_id"]
            initial_guess_ids.append(int(source_id))


    # sort by cluster! (see parent structure (reference id of nanocluster))
    # within reference_id_dct

    # get structures of parent nanoclusters #
    ########################################
    ref_ids = [ref_id for ref_id, value in reference_id_dct.items()]

    ref_simulations, ref_simulations_list = fetch_simulations(ref_ids, is_test = is_test)

    # get positions of initial adsorbate guess #
    ############################################
    initial_guess_simulations, initial_guess_simulations_list = fetch_simulations(initial_guess_ids, is_test = is_test)

    # All necessary simulations are downloaded #
    ############################################
    positions_dct = {}
    adstypes_dct = {}
    compositions = []
    composition_dct = {}
    y_property_shortened = []

    for latest_simulation, simulation in zip(simulations_list, initial_guess_simulations_list):
        operations = simulation["operations"]
        if operations != target_operations:
            print("WARNING, not add adsorbate operation:",  latest_simulation["_id"])
        # compositions
        composition = int(simulation["atoms"]["numbers"].count(element1))
        compositions.append(composition)

        value = composition_dct.get(composition, [])
        value.append(int(simulation["_id"]))
        composition_dct[composition] = value
        
        # atomic adsorbate positions, types
        atom_ids = simulation["adsorbates"][0]["atom_ids"]
        adspos = np.array(simulation["atoms"]["positions"])[np.array(atom_ids)]
        adstypes = np.array(simulation["atoms"]["numbers"])[np.array(atom_ids)]

        positions_dct[latest_simulation["_id"]] = adspos
        adstypes_dct[latest_simulation["_id"]] = adstypes

        # shortened list of adsorption energies (possibly the workflow contains more than two different elements)
        y_property_shortened.append(property_dct[latest_simulation["_id"]])


    compositions = np.array(compositions)

    # summarized structures of nanoclusters fully covered with adsorbate #
    ######################################################################

    covered_structures = []
    nanoclusters = []

    # sort nanoclusters by element 1
    nel1_lst = []
    for nanocluster in ref_simulations_list:
        nano_atoms_dct = nanocluster["atoms"]
        nano_atoms = atoms_dict_to_ase(nano_atoms_dct)
        nel1 = (nano_atoms.get_atomic_numbers() == element1).sum()
        nel1_lst.append(nel1)
    sort_ids = np.argsort(np.array(nel1_lst))
    ref_simulations_list = np.array(ref_simulations_list)[sort_ids]


    #print("ref_simulations_list", ref_simulations_list)
    for nanocluster in ref_simulations_list:
        ref_id = nanocluster["_id"]
        nano_atoms_dct = nanocluster["atoms"]
        nano_atoms = atoms_dict_to_ase(nano_atoms_dct)
        covered_atoms = nano_atoms.copy()

        related_ids = reference_id_dct[ref_id]
        for related_id in related_ids:
            positions = positions_dct[related_id]
            adstypes = adstypes_dct[related_id]
            y = property_dct[related_id]

            # using charges to store info
            adsatoms = ase.Atoms(numbers = adstypes, positions = positions,
                masses = [y] * len(adstypes), tags = [int(related_id)] * len(adstypes))
            covered_atoms += adsatoms

        nanoclusters.append(nano_atoms)
        covered_structures.append(covered_atoms.copy())

    # plot structures with adsorbate energies by color #
    ####################################################
    # get separate color bar mapping to adsorption energy (global)
    minimum = y_property.min()
    maximum = y_property.max()
    #minimum = -1.8
    #maximum = 0.0

    # prepare plotting structures
    plt.close('all')
    fig = plt.figure(figsize=(len(nanoclusters) * 10, 10))
    #fig = plt.figure()
    outer = gridspec.GridSpec(1, len(nanoclusters), wspace=0.2, hspace=0.1)


    for i, (nanocluster, structure) in enumerate(zip(nanoclusters, covered_structures)):
        nel1 = (nanocluster.get_atomic_numbers() == element1).sum()
        nel2 = (nanocluster.get_atomic_numbers() == element2).sum()

        filename_root = "zz_test_"  + str(i) + "_" + pte[element1] + str(nel1) + pte[element2] + str(nel2) + "_nc"
        filename1 = plot_atoms_properties(nanocluster, rotation = '', 
            filename = filename_root, minimum = minimum, maximum = maximum, is_render_skipped = is_render_skipped)
        filename2 = plot_atoms_properties(nanocluster, rotation = '180y', 
            filename = filename_root + "_rot", minimum = minimum, maximum = maximum, is_render_skipped = is_render_skipped)
        
        filename_root = "zz_test_"  + str(i) + "_" + pte[element1] + str(nel1) + pte[element2] + str(nel2) + "_cov"
        filename3 = plot_atoms_properties(structure, rotation = '', 
            filename = filename_root, minimum = minimum, maximum = maximum, is_render_skipped = is_render_skipped)
        filename4 = plot_atoms_properties(structure, rotation = '180y', 
            filename = filename_root + "_rot", minimum = minimum, maximum = maximum, is_render_skipped = is_render_skipped)

        filenames = [filename1, filename2, filename3, filename4,]

        # save adsorbate energies into extended xyz format
        ase.io.write(filename_root + ".xyz", structure)

        # plot nanoclusters + covered nanoclusters (2 sides) in subplots

        inner = gridspec.GridSpecFromSubplotSpec(2, 2, 
            subplot_spec=outer[i], wspace=0.0, hspace=0.0)

        indexing_matrix = np.arange(4).reshape((2,2))
        for j, filename in enumerate(filenames):
            j0,j1 = np.where(indexing_matrix == j)
            j0, j1 = j0[0], j1[0]
            print(j0, j1)
            ax = plt.Subplot(fig, inner[j0, j1])
            #t = ax.text(0.5,0.5, 'outer=%d, inner=%d' % (i,j))
            #t.set_ha('center')
            im = ax.imshow(load_image(filename), interpolation='none',
                origin='lower',
                #extent=[-2, 4, -3, 2], clip_on=True
                )
            ax.set_xticks([])
            ax.set_yticks([])
            if j == 0:
                ax.set_title(pte[element1] + str(nel1) + pte[element2] + str(nel2))
            fig.add_subplot(ax)

    #fig.show()
    plt.savefig("overview_structures.png", dpi = 300)
    #plt.show()


    #plt.close('all')

    ### plotting compositions vs adsorption energy ###
    ##################################################
    rootname = "BinTM_H"
    element = pte[element1] + pte[element2]
    y_property = np.array(y_property_shortened)

    ### PtFe NH3 data requires special shift!
    #Eref = -11.738098438387794
    #EH2 = -1.16195386047558 * 0.5 
    #offset = Eref - EH2 
    #offset = 0 
    #ev = 27.12
    #y_property = (y_property - offset ) * ev

    x_ratio = compositions / 55.0
    print("compositions", compositions.shape)
    x_mean = sorted(list((set(x_ratio))))
    y_mean = []
    for c in x_mean:
        y_mean.append(np.mean(y_property[c == x_ratio]))

    print(x_mean)
    print(y_mean)

    if len(x_mean) <=1:
        print("Not enough different compositions,", str(len(x_mean)), ". Exiting")
        exit(0)

    label = "ratio " + pte[element1] + " : (" + pte[element1] + " + " +  pte[element2] + ")"
    plt.close('all')

    plot_scatter_mean(x_ratio, y_property, 
        x_mean, y_mean,
        label = label, 
        name = rootname + "_" + pte[element1] + pte[element2] + "_mean_scatter.png")

    # modify scatterplot to be interactive
    if is_noninter == True:
        print("skipping interactive scatterplot")
    else:
        interactive_plot_scatter_mean(x_ratio, y_property, x_mean, y_mean, label = label, name = "noname_scatter.png")
