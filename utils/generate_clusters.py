import sys
import numpy as np
import unittest
from ase.build import molecule
from ase.visualize import view
from ase.cluster.icosahedron import Icosahedron
from ase.io import write
atoms = Icosahedron('Cu', noshells=3)

from cluskit import Cluster
import cluskit
import dscribe

#cluster = Cluster(atoms)
cluster_list = []
path = "ptcu_clusters/"

# i in wulff is number of atoms
# in octahedron number of atoms on the edge
# in ico number of shells
for n_shells in [55,147,309]:
    #    for shape in ["ico", "octa", "wulff"]:
    for shape in ["wulff"]:

        # here energies, surfaces are for wulff construction (average over Cu and Pt from crystallium database)
        scaffold = cluskit.build.get_scaffold(shape = shape, i = n_shells, latticeconstant = 3.7,
            energies = [(1.47 + 1.84)/2.0, (1.31 + 1.48) / 2.0, (1.56 + 1.68) / 2.0], 
            surfaces = [(1, 0, 0), (1, 1, 1), (1, 1, 0)])

        # descriptor needs to be set correctly
        scaffold.descriptor_setup = dscribe.descriptors.SOAP(
            species=[29,78],
            periodic=False,
            rcut=5.0,
            nmax=8,
            lmax=6,
            sparse=False,
            average=True
            )

        # build special clusters
        cluster_list.extend( scaffold.get_segregated(typeA = 29, typeB = 78, ntypeB = 13, n_clus = 1) )
        write(path + "pt13cux_" + shape + "_" + str(n_shells) + "_seg" + ".xyz", cluster_list[-1].ase_object)

        cluster_list.extend( scaffold.get_core_shell(typeA = 29, typeB = 78, ntypeB = 13, n_clus = 1) )
        write(path + "pt13cux_" + shape + "_" + str(n_shells) + "_coreshell" + ".xyz", cluster_list[-1].ase_object)

        cluster_list.extend( scaffold.get_ordered(typeA = 29, typeB = 78, ntypeB = 13, n_clus = 1) )
        write(path + "pt13cux_" + shape + "_" + str(n_shells) + "_ordered" + ".xyz", cluster_list[-1].ase_object)

        cluster_list.extend( scaffold.get_random(typeA = 29, typeB = 78, ntypeB = 13, n_clus = 1) )
        write(path  + "pt13cux_" + shape + "_" + str(n_shells) + "_random" + ".xyz", cluster_list[-1].ase_object)


view(cluster_list)
