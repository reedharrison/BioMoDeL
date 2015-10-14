###################################################################################################################################
# Function:		calcElasticEnergy
# Purpose:		Calculate energy stored in spring deformations from elastic network model according to some reference structure
# Arguments:	
###################################################################################################################################
def calcElasticEnergy(pdb_ref, ensemble, sel, cutoff, k=1):
	import prody as pd
	import numpy as np
	from scipy import spatial

	# Select desired atoms, retrieve reference parameters, and initialize Kirchhoff matrix
	bb_ref = pdb_ref.select(sel)
	coords_ref = bb_ref.getCoords()
	gnm = pd.GNM()
	gnm.buildKirchhoff(bb_ref, cutoff, k)
	kirchhoff = gnm.getKirchhoff()
	np.fill_diagonal(kirchhoff, 0)
	kirchhoff = abs(kirchhoff)
	kirchhoff = spatial.distance.squareform(kirchhoff)

	# Calculate distances between Kirchhoff contacts
	eq_dists = spatial.distance.pdist(coords_ref)

	# Iterate over coords within ensemble
	energies = np.zeros(len(ensemble))
	ind = 0
	for coords in ensemble.iterCoordsets():
		print('Calculating elastic energies for element %d\n' % ind)
		dists = spatial.distance.pdist(coords)
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
	from scipy import spatial
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
		kirchhoff = spatial.distance.squareform(kirchhoff)

		# Calculate deformation energies
		eq_dists = spatial.distance.pdist(coords_ref)
		energies = np.zeros(len(ensemble))
		j = 0
		for coords in ensemble.iterCoordsets():
			dists = spatial.distance.pdist(coords)
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
# Purpose:		Generate coordinates for landscape structures. Each structure is linked to projections on PC axes for plotting.
# Arguments:	
###################################################################################################################################
def makeLandscape(pdb, ensemble, sel='name CA', nstep = 100, prefix = 'landscape'):
	import prody as pd
	import scipy as sci
	import numpy as np

	# Initialize ensemble
	atoms = pdb.select(sel)
	ensemble.setAtoms(atoms)
	num_atom = len(ensemble.getAtoms())

	# Calculate PCA modes of transition path
	pca = pd.PCA('PCA of Transition Path')
	pca.buildCovariance(ensemble)
	pca.calcModes()

	# Generate mesh for PC0 and PC1
	proj0 = pd.calcProjection(ensemble, pca[0])
	proj1 = pd.calcProjection(ensemble, pca[1])
	x = np.linspace(min(proj0),max(proj0),nstep)
	y = np.linspace(min(proj1),max(proj1),nstep)
	xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
	xv = xv.reshape((xv.size,1))
	yv = yv.reshape((yv.size,1))
	proj_land = np.column_stack((xv,yv))

	# Generate landscape structures
	coords_land = z2xyz(proj_land, ensemble, pca[(0,1)])
	coords_land = coords_land.reshape((coords_land.shape[0],coords_land.shape[1]/3,3))
	ensemble_land = pd.Ensemble()
	ensemble_land.addCoordset(coords_land)

	return((proj_land,ensemble_land))

###################################################################################################################################
# Function:		xyz2z
# Purpose:		Convert a z-score(s) from a principle component to the associated xyz vector. Number of principle components used
#				for xyz construction must match the number of eigenvectors provided
# Arguments:	
###################################################################################################################################
# def xyz2z(ensemble, modes, *args, **kwargs):
# 	proj = calcProjection(ensemble, modes, kwargs.pop('rmsd', True))


###################################################################################################################################
# Function:		z2xyz
# Purpose:		Convert a z-score(s) from a principle component to the associated xyz vector. Each column of z corresponds to a PC,
#				and the modes object must only contain information from that number of PCs
# Arguments:	
###################################################################################################################################
def z2xyz(z, ensemble, modes, rmsd = True):
	import numpy as np
	import prody as pd

	num_atom = ensemble.numAtoms()
	# deviations = ensemble.getDeviations()
	# if deviations.ndim == 3:
	# 	deviations = deviations.reshape((deviations.shape[0], deviations.shape[1] * 3))
	# elif deviations.ndim == 2:
	# 	deviations = deviations.reshape((1, deviations.shape[0] * 3))
	# else:
	# 	deviations = deviations.reshape((1, deviations.shape[0]))
	if rmsd:
		z = (num_atom**0.5) * z
	evec = modes.getArray()

	# Make sure shape of evec is properly recognized
	try:
		evec.shape[0]
	except:
		shape0 = 1
	else:
		shape0 = evec.shape[0]
	try:
		evec.shape[1]
	except:
		shape1 = 1
	else:
		shape1 = evec.shape[1]
	evec = evec.reshape((shape0, shape1))

	# Make sure shape of z is properly recognized
	try:
		z.shape[0]
	except:
		shape0 = 1
	else:
		shape0 = z.shape[0]
	try:
		z.shape[1]
	except:
		shape1 = 1
	else:
		shape1 = z.shape[1]
	z = z.reshape((shape0, shape1))

	# Calculate coordinates from supplied PCs
	xyz = ensemble[0].getCoords()
	xyz = xyz.reshape((1,xyz.shape[0]*xyz.shape[1]))
	coords = np.dot(evec, z.transpose())

	coords = coords.transpose() + xyz

	return(coords)

###################################################################################################################################
# Function:		f_hookean
# Purpose:		Equation for calculating elastic energies
# Arguments:	eq_dists = equilibrium distances in ndarray, index must correspond to distance between node i and j
#				dists = non-equilibrium distances in ndarray
#				k = spring constant in kirchhoff matrix (0 if no spring for interaction i-j)
###################################################################################################################################
def f_hookean(eq_dists, dists, k):
	import numpy as np
	from scipy import spatial
	return(np.sum(spatial.distance.squareform(k*(dists-eq_dists)**2)))