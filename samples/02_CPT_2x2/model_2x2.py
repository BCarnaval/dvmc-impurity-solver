"""This module aims to test the link(s) between dVMC and PyQCM
by using the new interface of PyQCM for CPT calculation.
"""

import pyqcm


# Global attributes of the cluster
ns, nb = 4, 0
no = ns + nb
CM = pyqcm.cluster_model(name='2x2', n_sites=ns, n_bath=nb)

# Construction of the lattice model
cluster = pyqcm.cluster(X=CM, sites=[[0, 0, 0], [1, 0, 0],
                                     [0, 1, 0], [1, 1, 0]])
model = pyqcm.lattice_model(
    name='2x2', clus=cluster, superlattice=[[2, 0, 0], [0, 2, 0]])

model.interaction_operator('U')
model.hopping_operator('t', [1, 0, 0], -1)  # NN hopping
model.hopping_operator('t', [0, 1, 0], -1)  # NN hopping
model.hopping_operator('tp', [1, 1, 0], -1)  # NNN hopping
model.hopping_operator('tp', [-1, 1, 0], -1)  # NNN hopping
model.hopping_operator('tpp', [2, 0, 0], -1)  # NNN hopping
model.hopping_operator('tpp', [0, 2, 0], -1)  # NNN hopping
