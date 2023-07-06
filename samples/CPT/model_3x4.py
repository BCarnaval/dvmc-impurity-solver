"""This module aims to test the link(s) between dVMC and PyQCM
by using the new interface of PyQCM for CPT calculation.
"""

import pyqcm


ns, nb = 12, 0
no = ns + nb
CM = pyqcm.cluster_model(name='3x4', n_sites=ns, n_bath=nb)

# 2x2_8b, tb1, one-body operator
# CM.new_operator(name='tb1', op_type='one-body', elem=[
#    (1, 5, -1.0),
#    (2, 6, -1.0),
#    (3, 7, -1.0),
#    (4, 8, -1.0),
#    (1+no, 5+no, -1.0),
#    (2+no, 6+no, -1.0),
#    (3+no, 7+no, -1.0),
#    (4+no, 8+no, -1.0)
# ])

# 2x2_8b, tb2, one-body operator
# CM.new_operator(name='tb2', op_type='one-body', elem=[
#    (1, 9, -1.0),
#    (2, 10, -1.0),
#    (3, 11, -1.0),
#    (4, 12, -1.0),
#    (1+no, 9+no, -1.0),
#    (2+no, 10+no, -1.0),
#    (3+no, 11+no, -1.0),
#    (4+no, 12+no, -1.0)
# ])

# 2x2_8b, eb1, one-body operator
# CM.new_operator(name='eb1', op_type='one-body', elem=[
#    (5, 5, 1.0),
#    (6, 6, 1.0),
#    (7, 7, 1.0),
#    (8, 8, 1.0),
#    (5+no, 5+no, 1.0),
#    (6+no, 6+no, 1.0),
#    (7+no, 7+no, 1.0),
#    (8+no, 8+no, 1.0)
# ])

# 2x2_8b, eb2, one-body operator
# CM.new_operator(name='eb2', op_type='one-body', elem=[
#    (9, 9, 1.0),
#    (10, 10, 1.0),
#    (11, 11, 1.0),
#    (12, 12, 1.0),
#    (9+no, 9+no, 1.0),
#    (10+no, 10+no, 1.0),
#    (11+no, 11+no, 1.0),
#    (12+no, 12+no, 1.0)
# ])


# Construction of the lattice model
cluster = pyqcm.cluster(X=CM, sites=[[0, 0, 0], [1, 0, 0], [2, 0, 0],
                                     [0, 1, 0], [1, 1, 0], [2, 1, 0],
                                     [0, 2, 0], [1, 2, 0], [2, 2, 0],
                                     [0, 3, 0], [1, 3, 0], [2, 3, 0]])
model = pyqcm.lattice_model(
    name='3x4', clus=cluster, superlattice=[[3, 0, 0], [0, 4, 0]])

model.interaction_operator('U')
model.hopping_operator('t', [1, 0, 0], -1)  # NN hopping
model.hopping_operator('t', [0, 1, 0], -1)  # NN hopping
model.hopping_operator('tp', [1, 1, 0], -1)  # NNN hopping
model.hopping_operator('tp', [-1, 1, 0], -1)  # NNN hopping
model.hopping_operator('tpp', [2, 0, 0], -1)  # NNN hopping
model.hopping_operator('tpp', [0, 2, 0], -1)  # NNN hopping
