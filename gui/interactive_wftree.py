# simple graph visualizer of the workflow
# TODO: location of nodes can be improved
# interactive, view structure upon clicking

import networkx as nx
import matplotlib.pyplot as plt 
from grave import plot_network
from grave.style import use_attributes
from critcatworks.database.format import atoms_dict_to_ase, ase_to_atoms_dict
from critcatworks.database.extdb import fetch_simulations
import ase, ase.io
from ase.visualize import view
import pymongo
from pprint import pprint as pp
import argparse
import numpy as np
#import matplotlib.pyplot as plt 
#import matplotlib

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id')
    parser.add_argument('--test', action='store_true')

    args = parser.parse_args()
    print(args.ids)
    is_test = args.test
    for idx in args.ids:
        simulations = get_simulations_from_workflow(idx, is_test)
        dol, simulations = get_dol_from_simulations(simulations)
        print(dol)
        draw_workflow(dol, simulations)
