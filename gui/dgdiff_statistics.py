# coding: utf-8
import json
import numpy as np
from matplotlib import pyplot as plt
import os
import logging
from pprint import pprint as pp
import argparse
import pymongo 

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


def boltzmann_distribution(e, T=298, k=8.617333262145e-5):
    '''
    Computes the Boltzmann distribution for energies e.
    Param:
    e: numpy array of energies (default unit: eV)
    T: temperature in Kelvin
    k: Boltzmann constant (units should correspond to those of energies), default is eV/K
    '''
    kT_inv = 1.0/(k*T)
    # Compute the canonical partition function
    Q = 0
    for e_j in e:
        Q += np.exp(-e_j * kT_inv)
    p = np.exp(-e * kT_inv)/Q
    return p


#''' Boltzmann distribution test '''
#fig = plt.figure(figsize=(10,5))
#e = np.linspace(0, 0.5, 50) # E(structure i) - E(ground state), ARTIFICIAL VALUES
#for T in [298]:
#    p = boltzmann_distribution(e, T=T)
#    plt.ylabel('Boltzmann factor')
#    plt.xlabel('E(i) - E(ground state) [eV]')
#    plt.plot(e, p, 'or', label='T={}'.format(T))
#    plt.plot(e, p, 'k--')
#plt.title('Example of Boltzmann distribution for artificial energy values')
#plt.legend()
#plt.show()

# This method returns deltaG_diff_prob
def compute_deltaG_diff_prob(energies_i, energies_j, E_H2, T):
    '''
    energies_i: np.array of energies (eV), n adsorbed H atoms
    energies_j: np.array of energies (eV), n-1 adsorbed H atoms
    T: temperature for the Boltzmann distribution
    '''
    p_i_distr = boltzmann_distribution(energies_i-min(energies_i), T=T)
    p_j_distr = boltzmann_distribution(energies_j-min(energies_j), T=T)
    deltaG_diff_prob = 0
    normalization_factor = 0
    for e_i, p_i in zip(energies_i, p_i_distr):
        for e_j, p_j in zip(energies_j, p_j_distr):
            deltaG_ij = e_i - e_j - 0.5 * E_H2 + 0.24
            #print('deltaG_ij=Ei(n)-Ej(n-1)-0.5*E(H2)={}'.format(deltaG_ij))
            deltaG_ij_prob = deltaG_ij * p_i * p_j
            deltaG_diff_prob += deltaG_ij_prob
            #print('deltaG_ij_prob=deltaG_ij*pi*pj={}'.format(deltaG_ij_prob))
            normalization_factor += p_i * p_j
    #print('normalization_factor={}'.format(normalization_factor))
    #print('deltaG_diff_prob={} (before normalization)'.format(deltaG_diff_prob))
    deltaG_diff_prob = deltaG_diff_prob / normalization_factor
    #print('deltaG_diff_prob={} (after normalization)'.format(deltaG_diff_prob))
    return deltaG_diff_prob


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('ids', metavar='N', type=int, nargs='+',
                    help='workflow id')
    parser.add_argument('--test', action='store_true')

    args = parser.parse_args()
    print(args.ids)
    workflow_id = args.ids[0]
    is_test = args.test

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
pp(dct)
dct = dct["output"]


branch_dct = {}
ne_dct = {}
reference_energy = 0.0
free_energy_correction = 0.0
json_files = ["FW_wf_57.json"]
#json_files = ["FW_wf_47.json", "FW_wf_56.json"]
#with open(json_file) as f:
#    dct = json.load(f)

step_history = dct["step_history"]
root_history = dct["root_history"]
branch_dct = dct["branch_dct"]
ne_dct = dct["ne_dct"]
        
reference_energy = dct["reference_energy"]
free_energy_correction = dct["free_energy_correction"]

prev_k = 10000000
for k, v in ne_dct.items():
    if int(k) < prev_k:
        prev_k = int(k)
lowest_k = prev_k
n_adsorbates = []
mean_n = []
results = []
mean_dg = []
median_dg = []
for k, v in ne_dct.items():
    if (int(k) - lowest_k ) >= 1:
        dgdiff = np.array(list(v.values())).reshape(1,-1) - np.array(list(ne_dct[str(int(k) - 1)].values())).reshape(-1,1) - reference_energy + free_energy_correction
        dgdiff = dgdiff * 27.2114
        print(dgdiff.shape)
        n_adsorbates.append(np.ones(dgdiff.shape) * int(k))
        results.append(dgdiff)

        mean_n.append(int(k))
        mean_dg.append(dgdiff.mean())
        median_dg.append(np.median(dgdiff))
flattened_results = np.concatenate(results, axis=None)
flattened_n = np.concatenate(n_adsorbates, axis=None)
#plt.scatter(flattened_n, flattened_results, facecolors='none', edgecolors='blue')
#plt.scatter(mean_n, mean_dg, s=80, facecolors='r', edgecolors='r')
#plt.scatter(mean_n, median_dg, s=80, facecolors='g', edgecolors='g')
#plt.xlabel('# adsorbates')
#plt.ylabel(r'$\Delta \ G \ [eV]$')
#plt.tight_layout()
#plt.savefig("dgdiff_stats_47_56_coverage_ladder.png", dpi = 300)
#plt.savefig("dgdiff_stats_57_coverage_ladder.png", dpi = 300)

#plt.show()
#plt.close()

# This is how I plot the figure
fig = plt.figure(figsize=(10,5))
deltaG_arr = []
T = 300
mean_n = []
ground_state_dgdiff = {}
ground_state_dgdiff_lst = []
for k, v in ne_dct.items():
    if (int(k) - lowest_k ) >= 1:
        energies_i = np.array(list(v.values())) * 27.2114
        energies_j = np.array(list(ne_dct[str(int(k) - 1)].values())) * 27.2114
        # NOTE: these are ALREADY in eV! and so is E_H2
        deltaG = compute_deltaG_diff_prob(energies_i, energies_j, reference_energy * 27.2114 * 2, T = 298)
        print('nH={}, deltaG_diff_prob={}'.format(mean_n,deltaG))
        deltaG_arr.append(deltaG)
        mean_n.append(int(k))
        # ground states
        dGdiff_ground = np.min(energies_i) - np.min(energies_j) - reference_energy * 27.2114 + free_energy_correction 
        ground_state_dgdiff[int(k)] = dGdiff_ground
        ground_state_dgdiff_lst.append(dGdiff_ground)
y_err = np.zeros(len(deltaG_arr))
# Plot with errorbars (todo...)
new_x, new_y, new_err = zip(*sorted(zip(mean_n, deltaG_arr, y_err)))
plt.errorbar(new_x, new_y, yerr=new_err, 
             color='r', marker='o', linestyle='--', markersize=5.0,
             ecolor='r', elinewidth=1, capsize=1, capthick=1,
             label='$\Delta G_\mathrm{diff}$ (all pairs + Boltzmann, 300 K)',
             uplims=True, lolims=True)
# Also plot ground state based deltaG values
y = np.array(ground_state_dgdiff_lst)
#plt.scatter(mean_n, y, c='b')
new_x, new_y = zip(*sorted(zip(mean_n, y)))
plt.plot(new_x, new_y, c='b', marker='o', markersize=5.0, linestyle='--', label='$\Delta G_\mathrm{diff}$ (ground states only)')
plt.plot([np.min(mean_n), np.max(mean_n)], [0, 0], 'k')
plt.xlabel('$n_\mathrm{H}$')
plt.ylabel('$\Delta G_\mathrm{diff}$ [eV]')
#plt.ylim([-0.2, 0.2])
#plt.xlim([np.min(mean_n), np.max(mean_n)])
plt.grid('on')
plt.legend()
plt.show()
filename = 'pt12ni42_deltaG_diff_prob.png'
filepath = os.path.join(os.getcwd(), filename)
fig.savefig(filepath, bbox_inches='tight')

