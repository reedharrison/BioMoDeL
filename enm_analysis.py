###################################################################################################################################
# Function:		calcElasticEnergy
# Purpose:		Calculate energy stored in spring deformations from elastic network model according to some reference structure
# Arguments:	
###################################################################################################################################
def calcElasticEnergy(pdb, refpdb, sel='calpha', cutoff=12, k=1):
	import prody as pd
	import numpy as np
	import scipy as sci
	import scipy.spatial as sp

	# Select desired atoms, retrieve reference parameters, and initialize Kirchhoff matrix
	gnm = pd.GNM()
	gnm.buildKirchhoff(refpdb.calpha, cutoff, k)
	kirchhoff = gnm.getKirchhoff()
	np.fill_diagonal(kirchhoff, 0)
	kirchhoff = abs(kirchhoff)

	# Calculate distances between Kirchhoff contacts
	eq_dists = sp.distance.cdist(refpdb.calpha.select(sel).getCoords(), refpdb.calpha.select(sel).getCoords(), metric='euclidean')
	dists = sp.cdist(pdb.calpha.select(sel).getCoords(), pdb.calpha.select(sel).getCoords(), metric='euclidean')

	# Iterate over coords within ensemble
	energies = k * np.multiply(np.square(dists-eq_dists), kirchhoff)
	energies = sci.triu(energies)
	energies = np.sum(energies)
	
	return(energies)

###################################################################################################################################
# Function:		calcElasticEnergy2
# Purpose:		Calculate energy stored in spring deformations from elastic network model according to some reference structure
# Arguments:	
###################################################################################################################################
def calcElasticEnergy2(pdb_ref, ensemble, sel, cutoff, k=1):
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
def calcMultiStateEnergy(ensemble, ensemble_ref, cutoff=12, k=1, U0=None):
	import prody as pd
	import scipy as sci
	import scipy.spatial as sp
	import numpy as np

	n_conf = ensemble.numConfs()
	n_ref = ensemble_ref.numConfs()
	energies = np.zeros((n_ref, n_conf))

	for i, conf_ref in zip(xrange(n_ref), ensemble_ref.iterCoordsets()):
		gnm = pd.GNM()
		gnm.buildKirchhoff(conf_ref, cutoff, k)
		kirchhoff = gnm.getKirchhoff()
		np.fill_diagonal(kirchhoff, 0)
		kirchhoff = abs(kirchhoff)
		eq_dists = sp.distance.squareform(sp.distance.pdist(conf_ref, metric='euclidean'))
		for j, conf in zip(xrange(n_conf), ensemble.iterCoordsets()):
			dists = sp.distance.squareform(sp.distance.pdist(conf, metric='euclidean'))
			springs = (k/2) * np.multiply(np.square(dists-eq_dists), kirchhoff)
			springs = sci.triu(springs)
			energies[i,j] = np.sum(springs)

	# return np.min(energies, axis=0)
	return energies
	


#def calcMultiStateEnergy(pdb_sel, ensemble_ref, ensemble, sel, cutoff, k=1, U_0=0):
	# # Select specified atoms for analysis
	# atoms = pdb_sel.select(sel)
	# ensemble_ref.setAtoms(atoms)
	# ensemble.setAtoms(atoms)
	# num_ref = len(ensemble_ref)
	# num_struct = len(ensemble)

	# # Initialize containers
	# U = np.zeros((num_ref,num_struct))
	
	# # If user defines a U_0, use the values. U_0[i] is the reference state energy for reference structure i
	# if U_0 == 0:
	# 	U_0 = np.zeros(num_ref)

	# # Perform elastic energy calculations for each reference structure
	# i = 0
	# for coords_ref in ensemble_ref.iterCoordsets():
	# 	# Define connectivity matrix for reference structure i
	# 	gnm = pd.GNM()
	# 	gnm.buildKirchhoff(coords_ref, cutoff, k)
	# 	kirchhoff = gnm.getKirchhoff()
	# 	np.fill_diagonal(kirchhoff, 0)
	# 	kirchhoff = abs(kirchhoff)
	# 	kirchhoff = spatial.distance.squareform(kirchhoff)

	# 	# Calculate deformation energies
	# 	eq_dists = spatial.distance.pdist(coords_ref)
	# 	energies = np.zeros(len(ensemble))
	# 	j = 0
	# 	for coords in ensemble.iterCoordsets():
	# 		dists = spatial.distance.pdist(coords)
	# 		U[i,j] = f_hookean(eq_dists, dists, k*kirchhoff)
	# 		j += 1
	# 	U = U + U_0[i]
	# 	i += 1
	# return(U.min(axis=0)) # Return minimum energy per structure by column


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

###################################################################################################################################
# Function:		plot_landscape
# Purpose:		Plot energy landscape as a surface
# Arguments:	
###################################################################################################################################
def plot_landscape(pc1, pc2, energy, label=['PC1', 'PC2', 'Elastic Energy'], rstride=1, cstride=1):
	import matplotlib.pyplot as plt
	from matplotlib import cm
	from mpl_toolkits.mplot3d import Axes3D

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(pc1, pc2, energy, rstride=rstride, cstride=cstride, cmap=cm.coolwarm, linewidth=0)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	ax.set_xlabel(label[0])
	ax.set_ylabel(label[1])
	ax.set_zlabel(label[2])
	return((fig, ax, surf))

###################################################################################################################################
# Function:		calcANMPathway
# Purpose:		Calculate ensemble of structures predicting conformational transition between pdb_a and pdb_b
# Arguments:	
###################################################################################################################################
def calcANMPathway(pdb_a, pdb_b, k=0.1, r_c=15, U0_a=0, U0_b=0, sa=0.8, sb=0.4, t_rmsd=0.1, tol=10**(-4), crit_rmsd=1, m=100, max_iter=100):
	import numpy as np
	import prody as pd
	import scipy as sci
	import scipy.spatial as sp

	def calc_dU(coords, coords_ref, cutoff=r_c, k=k):
		gnm = pd.GNM()
		gnm.buildKirchhoff(coords_ref, cutoff, k)
		kirchhoff = gnm.getKirchhoff()
		np.fill_diagonal(kirchhoff, 0)
		kirchhoff = abs(kirchhoff)

		n_atom = coords.shape[0]

		xi = coords[:,0]
		yi = coords[:,1]
		zi = coords[:,2]
		xj = coords_ref[:,0]
		yj = coords_ref[:,1]
		zj = coords_ref[:,2]

		xi, xj = np.meshgrid(xi, xj)
		yi, yj = np.meshgrid(yi, yj)
		zi, zj = np.meshgrid(zi, zj)
		mag = np.sqrt(np.square(xi-xj)+np.square(yi-yj)+np.square(zi-zj))
		np.fill_diagonal(mag, -1)	

		D = sp.distance.squareform(sp.distance.pdist(coords, metric='euclidean'))
		D0 = sp.distance.squareform(sp.distance.pdist(coords_ref, metric='euclidean'))

		dU = np.multiply(kirchhoff, D-D0)
		dU = dU/np.max(abs(dU))
		dUx = np.multiply(dU, np.divide(xi-xj, mag))
		dUy = np.multiply(dU, np.divide(yi-yj, mag))
		dUz = np.multiply(dU, np.divide(zi-zj, mag))

		# dUx = np.nansum(sci.triu(dUx))
		# dUy = np.nansum(sci.triu(dUy))
		# dUz = np.nansum(sci.triu(dUz))
		dUx = np.sum(dUx, axis=1)/dU.shape[1]
		dUy = np.sum(dUy, axis=1)/dU.shape[1]
		dUz = np.sum(dUz, axis=1)/dU.shape[1]

		return dUx, dUy, dUz

	def findCuspStruct(pdb_a, pdb_b, ensemble_ref, m=m):
		ensemble = pd.Ensemble()
		ensemble.setAtoms(pdb_a)
		ensemble.setCoords(pdb_a.getCoords())
		conf_i = pdb_a.copy()
		conf_f = pdb_b.copy()
		conf_f, T = pd.superpose(conf_f, conf_i)
		v = conf_f.getCoords() - conf_i.getCoords()
		for i in np.linspace(0, 1, m):
			q = i
			p = 1-q
			coords = (p*v)+conf_i.getCoords()
			ensemble.addCoordset(coords)
		E_trans = calcMultiStateEnergy(ensemble, ensemble_ref, cutoff=r_c, k=k)
		E_trans = E_trans/np.max(E_trans)
		diff_E = abs(E_trans[0,:]-E_trans[1,:])
		ind_trans = np.argmin(diff_E)
		coords = ensemble[ind_trans].getCoords()
		return(coords, diff_E[ind_trans])

	def minimize(coords, coords_ref, s, cutoff=r_c, k=k, U0=None):
		dUx, dUy, dUz = calc_dU(coords, coords_ref, cutoff=cutoff, k=k)
		dx = np.multiply(s, dUx)
		dy = np.multiply(s, dUy)
		dz = np.multiply(s, dUz)
		# print '\tMoving coordinates max <%f, %f, %f>'%(np.max(abs(dx)),np.max(abs(dy)),np.max(abs(dz)))
		x = coords[:,0] - dx
		y = coords[:,1] - dy
		z = coords[:,2] - dz
		newcoords = np.zeros(coords.shape)
		newcoords[:,0] = x
		newcoords[:,1] = y
		newcoords[:,2] = z
		return newcoords

	# Instantiate containers for data
	# pdb_b, junk = pd.superpose(pdb_b, pdb_a)
	pdb_container_a = pdb_a.copy()
	pdb_container_b = pdb_b.copy()
	pdb_trans = pdb_a.copy()
	path_a = pd.Ensemble('Path from transition to state A')
	path_b = pd.Ensemble('Path from transition to state B')
	path = pd.Ensemble('Transition Path')
	path_a.setAtoms(pdb_a)
	path_b.setAtoms(pdb_b)
	path.setAtoms(pdb_trans)
	# path_a.setCoords(pdb_a)
	# path_b.setCoords(pdb_b)

	ensemble_ref = pd.Ensemble()
	ensemble_ref.setAtoms(pdb_a)
	ensemble_ref.addCoordset(pdb_a)
	ensemble_ref.addCoordset(pdb_b)

	# Interpolate coordinates
	print 'Searching for initial transition state.'
	coords_trans_i, E_trans_i = findCuspStruct(pdb_container_a, pdb_container_b, ensemble_ref)

	# Search for transition state
	print 'Minimizing transition state.'
	coords_trans_f = coords_trans_i
	E_trans_f = E_trans_i
	counter = np.zeros(1)
	while((counter<max_iter)and(E_trans_f>tol)):
		counter += 1
		coords_trans_a = minimize(coords_trans_f, pdb_a.getCoords(), s=sa)
		coords_trans_b = minimize(coords_trans_f, pdb_b.getCoords(), s=sb)
		pdb_container_a.setCoords(coords_trans_a)
		pdb_container_b.setCoords(coords_trans_b)
		coords_trans_f, E_trans_f = findCuspStruct(pdb_container_a, pdb_container_b, ensemble_ref)
		print '\tBeginning iteration %d, dE=%f'%(counter, E_trans_f)
	pdb_trans.setCoords(coords_trans_f)

	# Find path from transition state to reference state A, using steepest descent
	print 'Finding paths of steepest descent from transition state.'
	counter = np.zeros(1)
	rmsd = pd.calcRMSD(pdb_a.getCoords(), pdb_trans.getCoords())
	pdb_container_a.setCoords(pdb_trans)
	while((counter<max_iter)and(rmsd>crit_rmsd)):
		counter += 1
		path_a.addCoordset(minimize(pdb_container_a.getCoords(), pdb_a.getCoords(), s=sa))
		pdb_container_a.setCoords(path_a[-1])
		rmsd = pd.calcRMSD(pdb_a.getCoords(), pdb_container_a.getCoords())
		print 'RMSD (path A): %f'%(rmsd)

	# Find path from transition state to reference state B, using steepest descent
	counter = np.zeros(1)
	rmsd = pd.calcRMSD(pdb_b.getCoords(), pdb_trans.getCoords())
	pdb_container_b.setCoords(pdb_trans)
	while((counter<max_iter)and(rmsd>crit_rmsd)):
		counter += 1
		path_b.addCoordset(minimize(pdb_container_b.getCoords(), pdb_b.getCoords(), s=sb))
		pdb_container_b.setCoords(path_b[-1])
		rmsd = pd.calcRMSD(pdb_b.getCoords(), pdb_container_b.getCoords())
		print 'RMSD (path B): %f'%(rmsd)

	# Stitch together frames of path in proper order
	for i in reversed(xrange(0, len(path_a))):
		path.addCoordset(path_a[i].getCoords())
	path.addCoordset(pdb_trans.getCoords())
	for i in xrange(0,len(path_b)):
		path.addCoordset(path_b[i].getCoords())
	print 'Transition path calculation complete!'
	
	return (path, pdb_trans)

###################################################################################################################################
# Function:		calcElasticPathway
# Purpose:		Calculate ensemble of structures predicting conformational transition between pdb_a and pdb_b
# Arguments:	
###################################################################################################################################
def calcElasticPathway(pdb_a, pdb_b, k=0.1, r_c=15, U0_a=0, U0_b=0, tol=10**(-4), crit_rmsd=1, max_iter=100):
	import numpy as np
	import prody as pd
	import scipy as sci
	import scipy.spatial as sp
	import scipy.optimize as so



	return 0