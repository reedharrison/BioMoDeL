###################################################################################################################################
# Function:		calcElasticEnergy
# Purpose:		Calculate energy stored in spring deformations from elastic network model according to some reference structure
# Arguments:	
###################################################################################################################################
def calcElasticEnergy(pdb_ref, ensemble, sel, cutoff, k=1):
	import prody as pd
	import numpy as np
	import scipy as sci

	# Select desired atoms, retrieve reference parameters, and initialize Kirchhoff matrix
	bb_ref = pdb_ref.select(sel)
	coords_ref = bb_ref.getCoords()
	gnm = pd.GNM()
	gnm.buildKirchhoff(bb_ref, cutoff, k)
	kirchhoff = gnm.getKirchhoff()
	np.fill_diagonal(kirchhoff, 0)
	kirchhoff = abs(kirchhoff)
	kirchhoff = sci.spatial.distance.squareform(kirchhoff)

	# Calculate distances between Kirchhoff contacts
	eq_dists = sci.spatial.distance.pdist(coords_ref)

	# Iterate over coords within ensemble
	energies = np.zeros(len(ensemble))
	ind = 0
	for coords in ensemble.iterCoordsets():
		print('Calculating elastic energies for element %d\n' % ind)
		dists = sci.spatial.distance.pdist(coords)
		energies[ind] = f_hookean(eq_dists, dists, k*kirchhoff)
		ind += 1
	
	return(energies)

###################################################################################################################################
# Function:		calcMinElasticLandscape
# Purpose:		Calculate minimum energy landscape according to Avisek et al, Plos Comp Bio, 2014
# Arguments:	
###################################################################################################################################
def calcMultiStateEnergy(pdb_sel, ensemble_ref, ensemble, sel, cutoff, k=1, U_0=0):
	import prody as pd
	import scipy as sci
	import numpy as np

	# Select specified atoms for analysis
	atoms = pdb_sel.select(sel)
	ensemble_ref.setAtoms(atoms)
	ensemble.setAtoms(atoms)
	num_ref = len(ensemble_ref)
	num_struct = len(ensemble)

	# Initialize containers
	U = np.zeros((num_ref,num_struct))
	
	# If user defines a U_0, use the values. U_0[i] is the reference state energy for reference structure i
	if U_0 == 0:
		U_0 = np.zeros(num_ref)

	# Perform elastic energy calculations for each reference structure
	i = 0
	for coords_ref in ensemble_ref.iterCoordsets():
		# Define connectivity matrix for reference structure i
		gnm = pd.GNM()
		gnm.buildKirchhoff(coords_ref, cutoff, k)
		kirchhoff = gnm.getKirchhoff()
		np.fill_diagonal(kirchhoff, 0)
		kirchhoff = abs(kirchhoff)
		kirchhoff = sci.spatial.distance.squareform(kirchhoff)

		# Calculate deformation energies
		eq_dists = sci.spatial.distance.pdist(coords_ref)
		energies = np.zeros(len(ensemble))
		j = 0
		for coords in ensemble.iterCoordsets():
			dists = sci.spatial.distance.pdist(coords)
			U[i,j] = f_hookean(eq_dists, dists, k*kirchhoff)
			j += 1
		U = U + U_0[i]
		i += 1
	return(U.min(axis=0)) # Return minimum energy per structure by column


###################################################################################################################################
# Function:		calcDistMat
# Purpose:		Calculate distance matrix given a matrix describing connectivity and a coordset
# Arguments:	k = Kirchhoff matrix describing node connectivity, zero describes no connection
#				coords = n x 3 matrix of xyz coordinates
###################################################################################################################################
# def calcDistMat(k, coords):
# 	import prody as pd
# 	import numpy as np
# 	import scipy as sci

# 	dists = sci.spatial.distance.pdist(coords)
# 	dists = sci.spatial.distance.squareform(dists)

# 	val = np.nditer(k, flags=['multi_index'])
# 	while not val.finished:
# 		if val[0] == 0:
# 			dists[val.multi_index[0], val.multi_index[1]] = 0

# 	return(dists)

###################################################################################################################################
# Function:		reduceSymMat
# Purpose:		Reduce symmetric matrix to unique elements
# Arguments:	m = matrix (ndarray from numpy) to be reduced
###################################################################################################################################
# def reduceSymMat(m, method = 0):
# 	import numpy as np

# 	dim = np.shape(m)
# 	for i in xrange(0,dim[0]):
# 		for j in xrange(0,dim[1]):
# 			if i<j:
# 				m[i,j] = 0

# 	return(m[m.nonzero()])

###################################################################################################################################
# Function:		getPeriodicCell
# Purpose:		Calculate energy stored in spring deformations from elastic network model according to some reference structure
# Arguments:	
###################################################################################################################################

###################################################################################################################################
# Function:		makeLandscape
# Purpose:		Calculate energy stored in spring deformations from elastic network model according to some reference structure
# Arguments:	
###################################################################################################################################

###################################################################################################################################
# Function:		z2xyz
# Purpose:		Convert a z-score(s) from a principle component to the associated xyz vector. Number of principle components used
#				for xyz construction must match the number of eigenvectors provided
# Arguments:	
###################################################################################################################################

###################################################################################################################################
# Function:		f_hookean
# Purpose:		Equation for calculating elastic energies
# Arguments:	eq_dists = equilibrium distances in ndarray, index must correspond to distance between node i and j
#				dists = non-equilibrium distances in ndarray
#				k = spring constant in kirchhoff matrix (0 if no spring for interaction i-j)
###################################################################################################################################
def f_hookean(eq_dists, dists, k):
	import numpy as np
	import scipy as sci
	return(np.sum(sci.spatial.distance.squareform(k*(dists-eq_dists)**2)))