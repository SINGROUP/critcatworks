# coding: utf-8
import json
import numpy as np
from matplotlib import pyplot as plt

step_history = []
root_history = []
branch_dct = {}
ne_dct = {}
reference_energy = 0.0
free_energy_correction = 0.0
json_files = ["../proximity/FW_wf_57.json"]
#json_files = ["FW_wf_47.json", "FW_wf_56.json"]
lst_json_files  =[["../similarity/FW_wf_47.json", "../similarity/FW_wf_56.json"],["../proximity/FW_wf_57.json"]]
#json_files = ["FW_wf_57.json", "FW_wf_47.json", "FW_wf_56.json"]
#with open(json_file) as f:
#    dct = json.load(f)
ne_dcts = []
for json_files in lst_json_files:
    ne_dct = {}
    for json_file in json_files:
        with open(json_file) as f:
            dct = json.load(f)

        dct = dct["spec"]["temp"]
        step_history.extend(dct["step_history"])
        root_history.extend(dct["root_history"])

        for k, v in  dct["branch_dct"].items():
            branch_dct.setdefault(k, []).extend(v)
        
        for k, v in  dct["ne_dct"].items():
            ne_dct.setdefault(k, {}).update(v)
        
        reference_energy = dct["reference_energy "]
        free_energy_correction = dct["free_energy_correction"]
    ne_dcts.append(ne_dct)


    
#ne_dct = dct["spec"]["temp"]
#ne_dct = dct["spec"]["temp"]["ne_dct"]
min_toten_dct = {}
min_id_dct = {}
for ne_dct, c in zip(ne_dcts, ["blue", "red"]):
    n_adsorbates = []
    root_x = []
    root_y = []
    mean_n = []
    results = []
    mean_dg = []
    median_dg = []
    for k, v in ne_dct.items():
        toten = np.array(list(v.values()))
        ids = np.array(list(v.keys()))
        if c == "red":
            min_toten = np.min(toten)
            min_toten_dct[str(k)] = min_toten
            min_id_dct[str(k)] = ids[np.argmin(toten)]
        else:
            min_toten = min_toten_dct.get(str(k), np.min(toten))
        toten = toten -min_toten
        toten = toten * 27.2114
        toten = toten[toten < 10]
        print(toten.shape)
        n_adsorbates.append(np.ones(toten.shape) * int(k))
        results.append(toten)
        mean_n.append(int(k))
        mean_dg.append(toten.mean())
        median_dg.append(np.median(toten))
        for idx, toten in v.items():
            if int(idx) in root_history:
                en = (toten - min_toten) * 27.2114
                if en < 10:
                    root_x.append(int(k))
                    root_y.append(en)
    flattened_results = np.concatenate(results, axis=None)
    flattened_n = np.concatenate(n_adsorbates, axis=None)
    if c == "blue":
        offset = -0.05
    else:
        offset = 0.05
    plt.scatter(flattened_n + offset, flattened_results, facecolors='none', edgecolors=c)
    plt.scatter(np.array(root_x) + offset, root_y, facecolors='green', edgecolors="green")

print(min_id_dct)
# compare with olli's datapoint , human workflow
E = -8769.184687580
en = E - min_toten_dct["80"]
plt.scatter(80, en, s=80, facecolors='black', edgecolors='black')

#plt.scatter(mean_n, mean_dg, s=80, facecolors='r', edgecolors='r')
#plt.scatter(mean_n, median_dg, s=80, facecolors='g', edgecolors='g')
plt.xlabel('# adsorbates')
plt.ylabel(r'total energy (w.r.t lowest energy) [eV]')
plt.legend(["similarity", "root", "proximity", "root", "human workflow"])
plt.tight_layout()
#plt.savefig("dgdiff_stats_47_56_coverage_ladder.png", dpi = 300)
plt.savefig("toten_stats.png", dpi = 300)

plt.show()
