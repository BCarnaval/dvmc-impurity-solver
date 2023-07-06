"""This module defines the lattice model using PyQCM objects.
"""
import pyqcm


# Global attributes of the cluster
ns, nb = 2, 4
no = ns + nb
CM = pyqcm.cluster_model(name='L2_4b', n_sites=ns, n_bath=nb)

# Defining interaction/hopping operators
CM.new_operator(name='tb1', op_type='one-body', elem=[
    (1, 3, -1.0),
    (2, 4, -1.0),
    (7, 9, -1.0),
    (8, 10, -1.0)
])

CM.new_operator(name='tb2', op_type='one-body', elem=[
    (1, 5, -1.0),
    (2, 6, -1.0),
    (7, 11, -1.0),
    (8, 12, -1.0)
])

CM.new_operator(name='eb1', op_type='one-body', elem=[
    (3, 3, 1.0),
    (4, 4, 1.0),
    (9, 9, 1.0),
    (10, 10, 1.0)
])

CM.new_operator(name='eb2', op_type='one-body', elem=[
    (5, 5, 1.0),
    (6, 6, 1.0),
    (11, 11, 1.0),
    (12, 12, 1.0)
])

CM.new_operator(name='sb1', op_type='anomalous', elem=[
    (1, 3+no, -1.0),
    (2, 4+no, -1.0),
    (3, 1+no, 1.0),
    (4, 2+no, 1.0)
])

CM.new_operator(name='sb2', op_type='anomalous', elem=[
    (1, 5+no, -1.0),
    (2, 6+no, -1.0),
    (5, 1+no, 1.0),
    (6, 2+no, 1.0)
])

CM.new_operator(name='pb1', op_type='anomalous', elem=[
    (3, 4+no, 1.0),
    (4, 3+no, -1.0),
    (5, 6+no, 1.0),
    (6, 5+no, -1.0)
])

# Construction of the lattice model
cluster = pyqcm.cluster(X=CM, sites=[[0, 0, 0], [1, 0, 0]])
model = pyqcm.lattice_model(
    name='1D_2_4b', clus=cluster, superlattice=[[2, 0, 0]])

model.interaction_operator('U')
model.hopping_operator('t', [1, 0, 0], -1)  # NN hopping

# NN hopping with imaginary amplitude
model.hopping_operator('ti', [1, 0, 0], -1, tau=2)
model.hopping_operator('tp', [2, 0, 0], -1)  # NNN hopping
model.hopping_operator('sf', [0, 0, 0], -1, tau=0,
                       sigma=1)  # on-site spin flip
model.hopping_operator('h', [0, 0, 0], -1, tau=0, sigma=3)  # on-site spin flip
model.anomalous_operator('D', [1, 0, 0], 1)  # NN singlet
# NN singlet with imaginary amplitude
model.anomalous_operator('Di', [1, 0, 0], 1j)
model.anomalous_operator('S', [0, 0, 0], 1)  # on-site singlet

# On-site singlet with imaginary amplitude
model.anomalous_operator('Si', [0, 0, 0], 1j)
model.anomalous_operator('dz', [1, 0, 0], 1, type='dz')  # NN triplet
model.anomalous_operator('dy', [1, 0, 0], 1, type='dy')  # NN triplet
model.anomalous_operator('dx', [1, 0, 0], 1, type='dx')  # NN triplet
model.density_wave('M', 'Z', [1, 0, 0])
model.density_wave('pT', 'dz', [1, 0, 0], amplitude=1, link=[1, 0, 0])
