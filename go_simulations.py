###################################################################################################################################
# Function:		readRestart
# Purpose:		Import data from Go restart file as a PDB object
# Arguments:	
###################################################################################################################################
def parseRestart(filein):
	import sys
	import os
	import time as t
	import numpy as np

	startTime = t.time()

	fw_header = (-7, 21, -3)
	fw_coor = (8, 8, 20, 20, 20)
	parse_header = make_parser(fw_header)
	parse_coor = make_parser(fw_coor)

	p = open(filein, 'r')

	time = np.zeros(1)
	mol = np.empty(0)
	serial = np.empty(0)
	coord = np.empty(0)

	for line in p:
		if " time = " in line:
			data = parse_header(line)
			time[0] = float(data[0])
		if " time = " not in line:
			data = parse_coor(line)
			mol = np.append(mol, int(data[0]))
			serial = np.append(serial, int(data[1]))
			coord = np.append(coord, map(float, data[2:5]))

	p.close()

	coord = coord.reshape((len(coord)/3,3))

	endTime = t.time()

	print 'Restart file read in %.2f seconds' %(endTime-startTime)

	return((time, mol, serial, coord))

###################################################################################################################################
# Function:		make_parser
# Purpose:		Generate a function to read a line of text into fields of a specified number of characters.
# Arguments:	fieldwidths = list of numbers telling how many characters to expect in a field. negative numbers are ignored.
#								widths must be given in expected order.
###################################################################################################################################
def make_parser(fieldwidths):
	try:
		from itertools import izip_longest  # added in Py 2.6
	except ImportError:
		from itertools import zip_longest as izip_longest  # name change in Py 3.x

	try:
		from itertools import accumulate  # added in Py 3.2
	except ImportError:
		def accumulate(iterable):
			'Return running totals (simplified version).'
			total = next(iterable)
			yield total
			for value in iterable:
				total += value
				yield total

	cuts = tuple(cut for cut in accumulate(abs(fw) for fw in fieldwidths))
	pads = tuple(fw < 0 for fw in fieldwidths) # bool values for padding fields
	flds = tuple(izip_longest(pads, (0,)+cuts, cuts))[:-1]  # ignore final one
	parse = lambda line: tuple(line[i:j] for pad, i, j in flds if not pad)
	# optional informational function attributes
	parse.size = sum(abs(fw) for fw in fieldwidths)
	parse.fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
		for fw in fieldwidths)
	return parse