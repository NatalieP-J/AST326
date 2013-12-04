import astropy.io.fits as ast
import numpy as np

def fitExtraction(fname,max_val = True,choice = None,save = False):
	'''Run this with just the file name as an argument and the function
		will return the spectrum containing the maximum value.

		If max_val set to False, it will return an array where each element
		is the set of values corresponding to a single spectrum

		If choice is set to an integer, the selected spectrum will be returned

		If save is set to True, the spectrum will be saved to a file.

	'''
	hdulist = ast.open('{0}.fit'.format(fname))
	spectra = hdulist[0].data
	if save == False:
		if max_val == True and choice == None:
			maxes = []
			for i in range(len(spectra)):
				maxes.append(max(spectra[i]))
			for j in range(len(maxes)):
				if maxes[j]==max(maxes):
					return spectra[j][::-1]
					print j
		if max_val == False and choice == None:
			return spectra
		if choice != None:
			return spectra[choice][::-1]
	if save == True:
		if max_val == True and choice == None:
			maxes = []
			for i in range(len(spectra)):
				maxes.append(max(spectra[i]))
			for j in range(len(maxes)):
				if maxes[j]==max(maxes):
					print j
					np.savetxt('{0}_spectra{1}.dat'.format(fname,j),spectra[j][::-1],delimiter='\n') 
		if max_val == False and choice == None:
			np.savetxt('{0}_spectra.dat'.format(fname),spectra,delimiter='\n')
		if choice != None:
			np.savetxt('{0}_spectra{1}.dat'.format(fname,choice),spectra[choice][::-1],delimiter='\n')

def returnrange(fname,column):
	h = ast.open('{0}.fit'.format(fname))
	h = h[0].data
	col = []
	for i in range(len(h)):
		col.append(h[i][column])
	mean = np.mean(col)
	std = np.std(col)
	print mean,std
	return col


