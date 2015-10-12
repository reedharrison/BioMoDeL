







#################################################################
# OLD FUNCTIONS TO PERFORM ITERATIVE ENM BELOW ... MUST UPDATE! #
#################################################################

def pdb2dcd( pdb, coord ):
	import numpy as np
	import prody as pd

	natoms = pdb.numAtoms()
	ncoord = np.size( dcd )[0]

	ensemble = pd.Ensemble()
	ensemble.setCoords(dcd[0,:].reshape(natoms,3))
	ensemble = pd.Ensemble()

	for i in xrange( 1, ncoord ):
		ensemble.addCoordset(dcd[i,:].reshape(natoms,3))

	return ensemble

def subset_ndarray( x, row, col):
    import numpy as np
    y = np.ndarray([np.size(col),np.size(row)],dtype=x.dtype)
    for i in xrange(0,np.size(row),1):
        for j in xrange(0,np.size(col),1):
            indR = row[i]
            indC = col[j]
            y[i,j] = x[indR,indC]
    return y

def norm_vec( v, nrow ):
	# This function accepts a vector to be interpreted as a matrix of nrow x 3 dimensions
	# This function typically operates on xyz type vector of coordinates
	# Next the function normalizes the vecto
	import numpy as np
	vi = v.reshape((nrow,3))
	vf = np.zeros((nrow,3))
	for i in xrange(0,nrow):
		mag = ( vi[i,0]**2 + vi[i,1]**2 + vi[i,2]**2 )**(0.5)
		vf[i,0] = vi[i,0] / mag
		vf[i,1] = vi[i,1] / mag
		vf[i,2] = vi[i,2] / mag
	vf = vf.reshape(nrow*3)
	return vf

def gauss_vec( v ):
	# This function accepts a matrix of xyz vectors (2D array)
	# Generates 1 random xyz vector based on each column in xyz matrix
	import numpy as np
	import random

	import numpy as np
	import random

	len_xyz = np.shape(v)[0]

	# Calculate random vector
	means = np.zeros((len_xyz))
	sdevs = np.zeros((len_xyz))
	gauss_vec = np.zeros((len_xyz))
	for i in xrange( 0, len_xyz ):
	    means[i] = np.average(v[i,:])
	    sdevs[i] = np.std(v[i,:])
	    gauss_vec[i] = random.gauss(means[i],sdevs[i])   
	return gauss_vec


# def gnm_figure():
# 	import numpy as np
# 	import matplotlib.pyplot as plt



# 	return 0

# def anm_figure():
# 	import numpy as np
# 	import matplotlib.pyplot as plt


# 	return 0

# def batch_gnm( pdbs, selstr, lambda, cutoff ):
# 	from prody import *
# 	from pylab import *

	
# 	return 0

# def batch_anm( pdbs, selstr, lambda, cutoff ):
# 	from prody import *
# 	from pylab import *


# 	return 0


def iANM_simple( pdb_obj, prefix, selstr='calpha', dr=1, nsteps=5, const=1, cutoff=10, selmode=0, path_length=0.5 ):
	import prody as pd
	import pylab as pb
	import numpy as np
	import random as random

	# This functions performs iterative ANM analysis, focusing on the motions associated with the
	#	selected mode relative to the slowest mode (randomly selects +/- component each iteration).

	pdb = pdb_obj.select( selstr )
	natoms = pdb.numAtoms()

	filename_pdb = '%s.pdb' % prefix
	filename_dcd = '%s.dcd' % prefix
	pd.writePDB( filename_pdb, pdb )

	dcd = np.zeros( ( natoms*3, 1 ) )
	dcd = pdb.getCoords().reshape( natoms*3 )

	ensemble = pd.Ensemble()
	ensemble.setCoords(dcd.reshape(natoms,3))

	for i in xrange( 1, nsteps+1 ):
		print 'Calculating coordinates at step %d\n' % i

		# Load coordinates from previous step
		pdb.setCoords( dcd.reshape((natoms,3)) )

		# ANM analysis
		anm = pd.ANM()
		anm.buildHessian( pdb, cutoff=cutoff, gamma=const )
		H = anm.getHessian()
		# N = (np.shape(H)[1] / 3)
		if np.size(selmode) > 1:
			N = max(selmode) + 1
		if np.size(selmode) == 1:
			N = selmode + 1
		anm.calcModes( n_modes=N,zeros=False,turbo=True )
		eVec = anm.getEigvecs()
		eVal = anm.getEigvals()

		# Normalize eigenvectors and generate gaussian if needed
		if np.size(selmode) > 1:
			rvec = np.zeros((np.shape(eVec)[0],np.size(selmode)))
			for i in xrange ( 0, np.size(eVec[1]) ):
				rvec[:,i] = norm_vec( eVec[:,selmode[i]], natoms )
			rvec = gauss_vec(rvec)
		if np.size(selmode) == 1:
			if selmode == 0:
				rvec = norm_vec( eVec, natoms )
			if selmode > 0:
				rvec = norm_vec( eVec[:,selmode], natoms )
		# if np.size(selmode) > 1:
		# 	for i in xrange (0, np.size(eVec[1]) ):
		# 		rvec = norm_vec( eVec[:,selmode[i]], natoms )
		# 	rvec = gauss_vec(rvec[:,selmode])
		# if np.size(selmode) == 1:
		# 	rvec = norm_vec( eVec[:,selmode], natoms )			

		rmag = np.zeros((natoms,3))
		fluct = pd.calcSqFlucts(anm[selmode])
		fluct = fluct / max(fluct)
		rmag[:,0] = fluct
		rmag[:,1] = fluct
		rmag[:,2] = fluct
		rmag = rmag.reshape(natoms*3)
		# rmag = ((-1)**random.randint(0,1)) * path_length * rmag # Commented out to look at positive eigenvector!
		rmag = path_length * rmag

		rmag = rmag.reshape(natoms*3)
		rvec = rvec.reshape(natoms*3)

		rvec = np.multiply(rvec,rmag)

		# Predict new coords
		dcd = rvec + pdb.getCoords().reshape( natoms*3 )
		ensemble.addCoordset(dcd.reshape(natoms,3))

		pd.writeDCD( filename_dcd, ensemble )

	return ( pdb, ensemble )

def iENM_stochastic( pdb_obj, prefix, selstr='calpha', dr=1, nsteps=1, const=1, cutoff=10, selmode=-1, path_length=0.5 ):
	import prody as pd
	import pylab as pb
	import numpy as np
	import random as random

	pdb = pdb_obj.select( selstr )
	natoms = pdb.numAtoms()

	dcd = np.zeros( ( natoms*3, 1 ) )
	dcd = pdb.getCoords().reshape( natoms*3 )

	ensemble = pd.Ensemble()
	ensemble.setCoords(dcd.reshape(natoms,3))

	filename_pdb = '%s.pdb' % prefix
	filename_dcd = '%s.dcd' % prefix
	pd.writePDB( filename_pdb, pdb )

	for i in xrange( 1, nsteps+1 ):
		print 'Calculating coordinates at step %d\n' % i

		# Load coordinates from previous step
		pdb.setCoords( dcd.reshape(natoms,3) )

		# GNM analysis
		gnm = pd.GNM()
		gnm.buildKirchhoff( pdb, cutoff=cutoff, gamma=const )
		K = gnm.getKirchhoff()
		N = np.shape(K)[1]
		gnm.calcModes( n_modes=N,zeros=False,turbo=True )
		if selmode == -1:
			fluct = pd.calcSqFlucts( gnm ) / max( pd.calcSqFlucts( gnm ) )
		if selmode != -1:
			fluct = pd.calcSqFlucts( gnm[selmode] ) / max( pd.calcSqFlucts( gnm[selmode] ) )

		# ANM analysis
		anm = pd.ANM()
		anm.buildHessian( pdb, cutoff=cutoff, gamma=const )
		H = anm.getHessian()
		N = (np.shape(H)[1] / 3)
		anm.calcModes( n_modes=N,zeros=False,turbo=True )
		eVec = anm.getEigvecs()
		eVal = anm.getEigvals()

		# Normalize vectors
		for i in xrange (0, np.size(eVec[1]) ):
			eVec[:,i] = norm_vec( eVec[:,i], natoms )

		# Randomly sample gaussian distribution of vectors
		gVec = gauss_vec(eVec)

		# Weight motion along gVec by gnm square fluctuations
		rvec = gVec
		rmag = np.zeros((natoms,3))
		rmag[:,0] = fluct
		rmag[:,1] = rmag[:,0]
		rmag[:,2] = rmag[:,0]
		rmag = rmag.reshape(natoms*3)
		rmag = ((-1)**np.random.random_integers(0,1,np.shape(rmag))) * path_length * rmag
		# rmag = ((-1)**random.randint(0,1)) * rmag
		rmag = rmag.reshape(natoms*3)
		rvec = rvec.reshape(natoms*3)
		rvec = np.multiply(rvec,rmag)

		# Predict new coords
		dcd = rvec.reshape(natoms*3) + pdb.getCoords().reshape(natoms*3)
		ensemble.addCoordset(dcd.reshape(natoms,3))
		pd.writeDCD( filename_dcd, ensemble )

	return ( pdb, ensemble )

def iENM_unconstrained( pdb_obj, prefix, selstr='calpha', dr=1, nsteps=1, const=1, cutoff=10, path_length=1 ):
	import prody as pd
	import pylab as pb
	import numpy as np
	import random as random

	pdb = pdb_obj.select( selstr )
	natoms = pdb.numAtoms()

	dcd = np.zeros( ( natoms*3, 1 ) )
	dcd = pdb.getCoords().reshape( natoms*3 )

	ensemble = pd.Ensemble()
	ensemble.setCoords(dcd.reshape(natoms,3))

	filename_pdb = '%s.pdb' % prefix
	filename_dcd = '%s.dcd' % prefix
	pd.writePDB( filename_pdb, pdb )

	for i in xrange( 1, nsteps+1 ):
		print 'Calculating coordinates at step %d\n' % i

		# Load coordinates from previous step
		pdb.setCoords( dcd.reshape(natoms,3) )

		# # GNM analysis
		# gnm = pd.GNM()
		# gnm.buildKirchhoff( pdb, cutoff=cutoff, gamma=const )
		# K = gnm.getKirchhoff()
		# N = np.shape(K)[1]
		# gnm.calcModes( n_modes=N,zeros=False,turbo=True )
		# fluct = pd.calcSqFlucts( gnm ) / max( pd.calcSqFlucts( gnm ) )

		# ANM analysis
		anm = pd.ANM()
		anm.buildHessian( pdb, cutoff=cutoff, gamma=const )
		H = anm.getHessian()
		N = (np.shape(H)[1] / 3)
		anm.calcModes( n_modes=N,zeros=False,turbo=True )
		eVec = anm.getEigvecs()
		eVal = anm.getEigvals()

		# Normalize vectors
		for i in xrange (0, np.size(eVec[1]) ):
			eVec[:,i] = norm_vec( eVec[:,i], natoms )

		# Randomly sample gaussian distribution of vectors
		gVec = gauss_vec(eVec)

		# Sample motion after given path_length 
		rvec = gVec
		rmag = np.ones((natoms,3))
		rmag = ((-1)**random.randint(0,1)) * path_length * rmag
		rmag = rmag.reshape(natoms*3)
		rvec = rvec.reshape(natoms*3)
		rvec = np.multiply(rvec,rmag)

		# Predict new coords
		dcd = rvec.reshape(natoms*3) + pdb.getCoords().reshape(natoms*3)
		ensemble.addCoordset(dcd.reshape(natoms,3))
		pd.writeDCD( filename_dcd, ensemble )

	return ( pdb, ensemble )

def iENM_deformationE( pdb_obj, ref_pdb_obj, prefix, selstr='calpha', dr=1, nsteps=1, const=1, cutoff=10, path_length=1 ):
	import prody as pd
	import pylab as pb
	import numpy as np
	import random as random

	pdb = pdb_obj.select( selstr )
	natoms = pdb.numAtoms()

	dcd = np.zeros( ( natoms*3, 1 ) )
	dcd = pdb.getCoords().reshape( natoms*3 )

	ensemble = pd.Ensemble()
	ensemble.setCoords(dcd.reshape(natoms,3))

	filename_pdb = '%s.pdb' % prefix
	filename_dcd = '%s.dcd' % prefix
	pd.writePDB( filename_pdb, pdb )

	ref = ref_pdb_obj.select( selstr )

	for i in xrange( 1, nsteps+1 ):
		print 'Calculating coordinates at step %d\n' % i

		# Load coordinates from previous step
		pdb.setCoords( dcd.reshape((natoms,3)) )

		# ANM analysis
		anm = pd.ANM()
		anm.buildHessian( pdb, cutoff=cutoff, gamma=const )
		H = anm.getHessian()
		# N = (np.shape(H)[1] / 3)
		if np.size(selmode) > 1:
			N = max(selmode) + 1
		if np.size(selmode) == 1:
			N = selmode + 1
		anm.calcModes( n_modes=N,zeros=False,turbo=True )
		eVec = anm.getEigvecs()
		eVal = anm.getEigvals()

		# Normalize eigenvectors and generate gaussian if needed
		if np.size(selmode) > 1:
			rvec = np.zeros((np.shape(eVec)[0],np.size(selmode)))
			for i in xrange ( 0, np.size(eVec[1]) ):
				rvec[:,i] = norm_vec( eVec[:,selmode[i]], natoms )
			rvec = gauss_vec(rvec)
		if np.size(selmode) == 1:
			if selmode == 0:
				rvec = norm_vec( eVec, natoms )
			if selmode > 0:
				rvec = norm_vec( eVec[:,selmode], natoms )
		# if np.size(selmode) > 1:
		# 	for i in xrange (0, np.size(eVec[1]) ):
		# 		rvec = norm_vec( eVec[:,selmode[i]], natoms )
		# 	rvec = gauss_vec(rvec[:,selmode])
		# if np.size(selmode) == 1:
		# 	rvec = norm_vec( eVec[:,selmode], natoms )			

		rmag = np.zeros((natoms,3))
		fluct = pd.calcSqFlucts(anm[selmode])
		fluct = fluct / max(fluct)
		rmag[:,0] = fluct
		rmag[:,1] = fluct
		rmag[:,2] = fluct
		rmag = rmag.reshape(natoms*3)
		# rmag = ((-1)**random.randint(0,1)) * path_length * rmag # Commented out to look at positive eigenvector!
		rmag = path_length * rmag

		rmag = rmag.reshape(natoms*3)
		rvec = rvec.reshape(natoms*3)

		rvec = np.multiply(rvec,rmag)

		# Predict new coords
		dcd = rvec + pdb.getCoords().reshape( natoms*3 )
		ensemble.addCoordset(dcd.reshape(natoms,3))

		pd.writeDCD( filename_dcd, ensemble )