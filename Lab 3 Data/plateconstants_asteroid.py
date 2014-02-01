import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
import urllib as url
import string as str
from scipy.optimize import leastsq

#SAMPLE OUTPUT
'''
a11 = [[ 1.00050002]], 
a12 = [[-0.00659824]], 
x0 = [[ 1019.64337372]], 
with reduced chi squared = [[ 1.43784932]]
a21 = [[ 0.00626632]], 
a22 = [[ 0.99983162]], 
y0 = [[ 1028.3480059]], 
with reduced chi squared = [[ 1.98449692]]

'''

ax = []
ay = []


errx = []
erry = []

day = [19,20,22,23,28]

for urania in range(5):
	print 'Jan{0}'.format(urania)
	skyname = '30Urania-S001-R001-C001-r-{0}.fts'.format(urania)
	darkname = 'Dark-S001-R001-C001-B2-{0}.fts'.format(urania)
	flatname = 'AutoFlat-Dusk-r-Bin2-001-{0}.fts'.format(urania)
	s = pf.open(skyname)
	sky_time = s[0].header['exptime']
	sky_date = s[0].header['date']
	sky_hour = s[0].header['time-obs']
	#print 'Raw exposure {0}s'.format(sky_time)
	#print 'Raw date {0}'.format(sky_date)
	#print 'Raw time {0}'.format(sky_hour)
	sky = pf.getdata(skyname)
	sky /= sky_time
	s = pf.open(darkname)
	dark_time = s[0].header['exptime']
	#print 'Dark exposure {0}s'.format(dark_time)
	dark = pf.getdata(darkname)
	dark /= dark_time
	s = pf.open(flatname)
	flat_time = s[0].header['exptime']
	#print 'Flat exposure {0}s'.format(flat_time)
	flat = pf.getdata(flatname)
	flat /= flat_time
	f = flat-dark
	f = f/np.median(f)
	correct = ((sky-dark)/f)[::-1]

	radian = np.pi/180

	#function that calls the database
	def unso(radeg,decdeg,fovam):
		str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1'
		str2 = '&-c.ra={0:4.6f}&-c.dec={1:4.6f}&-c.bm={2:4.7f}/{3:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)
		str = str1 + str2
		f = url.urlopen(str)
		s = f.read()
		f.close()
		s1 = s.splitlines()
		s1 = s1[45:-1]
		name = np.array([])
		rad = np.array([])
		ded = np.array([])
		rmag = np.array([])

		for k in s1:
			kw = k.split('\t')
			name = np.append(name,kw[0])
			rad = np.append(rad,float(kw[1]))
			ded = np.append(ded,float(kw[2]))
			if kw[12] != '     ':
				rmag = np.append(rmag,float(kw[12]))
			else:
				rmag = np.append(rmag,np.nan)
		return name,rad,ded,rmag

	galaxy = '30Urania-S001-R001-C001-r-{0}.fts'.format(urania)
	s1 = pf.open(galaxy)

	#read off the pointing infromation from the file
	ras = s1[0].header['ra']
	des = s1[0].header['dec']
	#convert ra from hours to degrees
	radeg = 15*(float(ras[0:2])+float(ras[3:5])/60. + float(ras[6:])/3600.)
	dsgn = np.sign(float(des[0:3]))
	#convert dec from deg:arcmin:arcsec to degrees
	dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.

	#calculate field of view

	f = 3420 #focal length of the telescope in mm
	p = 0.018 #pixel size in mm (account for 2x2 binning)

	#central pixel
	x0 = 1024
	y0 = 1024

	fov = 2048./(f/p)
	fovam = (fov/radian)*60 #arcminute field of view
	fov = 35

	#star locations from catalogue
	name,rad,ded,rmag = unso(radeg,dedeg,fovam)

	#centre point in radians
	ra0 = radeg*radian 
	dec0 = dedeg*radian 

	#star locations in radians
	rad *= radian
	ded *= radian

	#convert radians (on the sphere of the sky) to standard coordinates
	def rad2pix(p0,ra,dec):
		ra0,dec0 = p0
		X = -(np.cos(dec)*np.sin(ra-ra0))/((np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0))+(np.sin(dec)*np.sin(dec0)))
		Y = -((np.sin(dec0)*np.cos(dec)*np.cos(ra-ra0))-(np.cos(dec0)*np.sin(dec)))/((np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0))+(np.sin(dec)*np.sin(dec0)))
		return X,Y

	#standard coordinates of catalogues
	X,Y = rad2pix([ra0,dec0],rad,ded)

	#shift to pixels for ideal camera
	shiftedX = f*(X/p)+x0
	shiftedY = f*(Y/p)+y0

	#limit to stars brighter than 13
	w = np.where(rmag < 13.)[0]

	#import located centroids
	x = np.loadtxt('Jan{0}_Xcentroid.txt'.format(day[urania]))
	y = np.loadtxt('Jan{0}_Ycentroid.txt'.format(day[urania]))

	#limite catalogue results to those specified magnitude
	shiftX = shiftedX[w]
	shiftY = shiftedY[w]

	#match up catalouge and centroids via measures of separations
	separations = np.array([])

	#catalogue matched values
	Xvals = np.array([])
	Yvals = np.array([])

	#centroid matched values
	xvals = np.array([])
	yvals = np.array([])
	for m in range(len(x)):
		separation = np.array([])
		#temp list for each centroid
		databaseX = np.array([])
		databaseY = np.array([])
		dataX = np.array([])
		dataY = np.array([])
		for k in range(len(shiftX)):
			xdiff = abs(x[m]-shiftX[k])
			ydiff = abs(y[m]-shiftY[k])
			R = np.sqrt(xdiff**2+ydiff**2)
			#find matches with separation less than 15 pixels
			if R < 70:
				point = [R,shiftX[k],shiftY[k],x[m],y[m]]
				separation = np.append(separation,point[0])
				databaseX = np.append(databaseX,point[1])
				databaseY = np.append(databaseY,point[2])
				dataX = np.append(dataX,point[3])
				dataY = np.append(dataY,point[4])
		#if there is only one within that radius, use that one
		if len(separation) == 1:
			for l in range(len(separation)):
				if databaseX[l] not in Xvals: #check that the catalogue value hasn't already been matched
					separations = np.append(separations,separation[l])
					Xvals = np.append(Xvals,databaseX[l])
					Yvals = np.append(Yvals,databaseY[l])
					xvals = np.append(xvals,dataX[l])
					yvals = np.append(yvals,dataY[l])
		#if there are more, find the closest
		if len(separation) != 1:
			for l in range(len(separation)):
				if databaseX[l] not in Xvals: #check that the catalogue value hasn't already been matched
					if separation[l] == min(separation):
						separations = np.append(separations,separation[l])
						Xvals = np.append(Xvals,databaseX[l])
						Yvals = np.append(Yvals,databaseY[l])
						xvals = np.append(xvals,dataX[l])
						yvals = np.append(yvals,dataY[l])

	plt.figure()

	for i in range(len(Xvals)):
		if Xvals[i] not in shiftX:
			print 'Error'
		if Yvals[i] not in shiftY:
			print 'Error'

	for i in range(len(separations)):
		xdiff = abs(xvals[i]-Xvals[i])
		ydiff = abs(yvals[i]-Yvals[i])
		R = np.sqrt(xdiff**2+ydiff**2)
		if R == separations[i]:
			xtemp = [xvals[i],Xvals[i]]
			ytemp = [yvals[i],Yvals[i]]
			plt.plot(xtemp,ytemp,'m',linewidth = 3) #plot lines showing matched points
		if R != separations[i]:
			print 'Error'

	#overplot centroids with catalogue values
	plt.plot(x,y,'b.',label = 'Centroids')
	plt.plot(shiftX,shiftY,'r.',label ='Catalogue')
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	plt.title('Matched Catalogue and Centroid Values')
	plt.legend(loc='best')
	plt.ylim(0,2048)
	plt.xlim(0,2048)
	#plt.show()

	xdiff = xvals-Xvals
	ydiff = yvals-Yvals

	plt.figure()
	plt.plot(xvals,xdiff,'bd',label='x differences')
	plt.plot(yvals,ydiff,'gd',label='y differences')
	plt.xlabel('x or y [pixel]')
	plt.ylabel('Offset')
	plt.title('X and Y Offset Between Centroids and Catalogue Stars')
	plt.legend(loc='best')
	plt.xlim(0,2048)
	#plt.show()

	#telescope scale
	m = f/p

	#convert matched catalogue values to standard coordinates
	X_centroid = (Xvals-x0)/m
	Y_centroid = (Yvals-y0)/m

	#construct B matrix
	B = []
	for i in range(len(X_centroid)):
		B.append([m*X_centroid[i],m*Y_centroid[i],1])
	B = np.matrix(B)

	#centroid vector
	xvecs = []
	for k in range(len(xvals)):
		xvec = np.matrix([xvals[k],yvals[k],1])
		xvecs.append(xvec.T)

	xvecs = np.array(xvecs)

	xi = np.matrix(xvecs[:,0])
	yi = np.matrix(xvecs[:,1])

	#calculate x and y constants, as well as reduced chi squared for each
	x_consts = (np.linalg.inv(B.T*B))*B.T*xi
	M = xi - B*x_consts
	x_chi = M.T*M
	dof = len(Xvals)-3
	red_chi_x = x_chi/dof
	#print 'a11 = {0}, \na12 = {1}, \nx0 = {2}, \nwith reduced chi squared = {3}'.format(x_consts[0],x_consts[1],x_consts[2],red_chi_x)
	y_consts = (np.linalg.inv(B.T*B))*B.T*yi
	M = yi - B*y_consts
	y_chi = M.T*M
	dof = len(Yvals)-3
	red_chi_y = y_chi/dof
	#print 'a21 = {0}, \na22 = {1}, \ny0 = {2}, \nwith reduced chi squared = {3}'.format(y_consts[0],y_consts[1],y_consts[2],red_chi_y)

	T = np.matrix([[m*x_consts[0],m*x_consts[1],x_consts[2]],
		[m*y_consts[0],m*y_consts[1],y_consts[2]],
		[0,0,1]])

	#check the transformation by switching catalogue stars with new constants
	xi_check = B*x_consts
	yi_check = B*y_consts

	#plot the residuals
	plt.figure()
	plt.plot(xi,xi-xi_check,'b.',label='x')
	plt.plot(yi,yi-yi_check,'g.',label='y')
	plt.xlabel('x or y [pixel]')
	plt.ylabel('Residual')
	plt.xlim(0,2048)
	plt.title('Residuals')
	plt.legend(loc='best')
	#plt.show()

	xi = np.array(xi)
	yi = np.array(yi)
	xi_check = np.array(xi_check)
	yi_check = np.array(yi_check)

	xi2 = []
	yi2 = []
	xi_check2 = []
	yi_check2 = []
	Xvals2 = []
	Yvals2 = []
	xvecs2 = []
	outliersx = []
	outliersy = []

	for i in abs(xi-xi_check):
		if float(i)>10:
			outliersx.append(i)

	for i in abs(yi-yi_check):
		if float(i)>10:
			outliersy.append(i)

	for j in range(len(xi_check)):
		if abs(xi-xi_check)[j] in outliersx or abs(yi-yi_check)[j] in outliersy:
			pass
		else:
			xi2.append(xi[j])
			xi_check2.append(xi_check[j])
			Xvals2.append(Xvals[j])
			yi2.append(yi[j])
			yi_check2.append(yi_check[j])
			Yvals2.append(Yvals[j])
			xvecs2.append(xvecs[j])

	xi = np.array(xi2)
	yi = np.array(yi2)
	xi_check = np.array(xi_check2)
	yi_check = np.array(yi_check2)
	Xvals = np.array(Xvals2)
	Yvals = np.array(Yvals2)
	xvecs = np.array(xvecs2)

	plt.figure()
	plt.plot(xi,xi-xi_check,'b.',label='x')
	plt.plot(yi,yi-yi_check,'g.',label='y')
	plt.xlabel('x or y [pixel]')
	plt.ylabel('Residual')
	plt.xlim(0,2048)
	plt.title('Residuals with Outliers Removed')
	plt.legend(loc='best')
	#plt.show()

	xerror = np.mean(abs(xi-xi_check))
	yerror = np.mean(abs(yi-yi_check))

	print 'Error in x position is {0} pixels'.format(xerror)
	print 'Error in y position is {0} pixels'.format(yerror)


	#convert matched catalogue values to standard coordinates

	Xvecs = []

	for i in range(len(Xvals)):
		Xvec = np.matrix([Xvals[i],Yvals[i],0]).T
		Xvecs.append(Xvec)

	capitalXs = []
	for i in range(len(Xvals)):
		capitalX = np.array((np.linalg.inv(T)*Xvecs[i]).T)
		capitalXs.append(capitalX)

	X_centroid = []
	Y_centroid = []

	for i in range(len(Xvals)):
		X_centroid.append(capitalXs[i][0][0])
		Y_centroid.append(capitalXs[i][0][1])

	X_centroid = np.array(X_centroid)
	Y_centroid = np.array(Y_centroid)


	#construct B matrix
	B = []
	for i in range(len(X_centroid)):
		B.append([m*X_centroid[i],m*Y_centroid[i],1])
	B = np.matrix(B)

	xi = np.matrix(xvecs[:,0])
	yi = np.matrix(xvecs[:,1])

	#calculate x and y constants, as well as reduced chi squared for each
	x_consts = (np.linalg.inv(B.T*B))*B.T*xi
	M = xi - B*x_consts
	x_chi = M.T*M
	dof = len(Xvals)-3
	red_chi_x = x_chi/dof
	#print 'a11 = {0}, \na12 = {1}, \nx0 = {2}, \nwith reduced chi squared = {3}'.format(x_consts[0],x_consts[1],x_consts[2],red_chi_x)
	y_consts = (np.linalg.inv(B.T*B))*B.T*yi
	M = yi - B*y_consts
	y_chi = M.T*M
	dof = len(Yvals)-3
	red_chi_y = y_chi/dof
	#print 'a21 = {0}, \na22 = {1}, \ny0 = {2}, \nwith reduced chi squared = {3}'.format(y_consts[0],y_consts[1],y_consts[2],red_chi_y)

	T = np.matrix([[m*x_consts[0],m*x_consts[1],x_consts[2]],
		[m*y_consts[0],m*y_consts[1],y_consts[2]],
		[0,0,1]])

	#check the transformation by switching catalogue stars with new constants
	xi_check = B*x_consts
	yi_check = B*y_consts

	xi = np.array(xi)
	yi = np.array(yi)
	xi_check = np.array(xi_check)
	yi_check = np.array(yi_check)

	xi2 = []
	yi2 = []
	xi_check2 = []
	yi_check2 = []
	Xvals2 = []
	Yvals2 = []
	xvecs2 = []
	outliersx = []
	outliersy = []

	for i in abs(xi-xi_check):
		if float(i)>5:
			outliersx.append(i)

	for i in abs(yi-yi_check):
		if float(i)>5:
			outliersy.append(i)

	for j in range(len(xi_check)):
		if abs(xi-xi_check)[j] in outliersx or abs(yi-yi_check)[j] in outliersy:
			pass
		else:
			xi2.append(xi[j])
			xi_check2.append(xi_check[j])
			Xvals2.append(Xvals[j])
			yi2.append(yi[j])
			yi_check2.append(yi_check[j])
			Yvals2.append(Yvals[j])
			xvecs2.append(xvecs[j])

	xi = np.array(xi2)
	yi = np.array(yi2)
	xi_check = np.array(xi_check2)
	yi_check = np.array(yi_check2)
	Xvals = np.array(Xvals2)
	Yvals = np.array(Yvals2)
	xvecs = np.array(xvecs2)

	xerror = np.mean(abs(xi-xi_check))
	yerror = np.mean(abs(yi-yi_check))

	print 'Error in x position is {0} pixels'.format(xerror)
	print 'Error in y position is {0} pixels'.format(yerror)

	#plot the residuals
	plt.figure()
	plt.plot(xi,xi-xi_check,'b.',label='x')
	plt.plot(yi,yi-yi_check,'g.',label='y')
	plt.xlabel('x or y [pixel]')
	plt.ylabel('Residual')
	plt.xlim(0,2048)
	plt.title('Second Iteration Residuals')
	plt.legend(loc='best')
	#plt.show()

	#convert matched catalogue values to standard coordinates

	Xvecs = []

	for i in range(len(Xvals)):
		Xvec = np.matrix([Xvals[i],Yvals[i],0]).T
		Xvecs.append(Xvec)

	capitalXs = []
	for i in range(len(Xvals)):
		capitalX = np.array((np.linalg.inv(T)*Xvecs[i]).T)
		capitalXs.append(capitalX)

	X_centroid = []
	Y_centroid = []

	for i in range(len(Xvals)):
		X_centroid.append(capitalXs[i][0][0])
		Y_centroid.append(capitalXs[i][0][1])

	X_centroid = np.array(X_centroid)
	Y_centroid = np.array(Y_centroid)

	#print X_centroid

	#construct B matrix
	B = []
	for i in range(len(X_centroid)):
		B.append([m*X_centroid[i],m*Y_centroid[i],1])
	B = np.matrix(B)

	xi = np.matrix(xvecs[:,0])
	yi = np.matrix(xvecs[:,1])

	#calculate x and y constants, as well as reduced chi squared for each
	x_consts = (np.linalg.inv(B.T*B))*B.T*xi
	M = xi - B*x_consts
	x_chi = M.T*M
	dof = len(Xvals)-3
	red_chi_x = x_chi/dof
	#print 'a11 = {0}, \na12 = {1}, \nx0 = {2}, \nwith reduced chi squared = {3}'.format(x_consts[0],x_consts[1],x_consts[2],red_chi_x)
	y_consts = (np.linalg.inv(B.T*B))*B.T*yi
	M = yi - B*y_consts
	y_chi = M.T*M
	dof = len(Yvals)-3
	red_chi_y = y_chi/dof
	#print 'a21 = {0}, \na22 = {1}, \ny0 = {2}, \nwith reduced chi squared = {3}'.format(y_consts[0],y_consts[1],y_consts[2],red_chi_y)

	T = np.matrix([[m*x_consts[0],m*x_consts[1],x_consts[2]],
		[m*y_consts[0],m*y_consts[1],y_consts[2]],
		[0,0,1]])

	#check the transformation by switching catalogue stars with new constants
	xi_check = B*x_consts
	yi_check = B*y_consts

	xi = np.array(xi)
	yi = np.array(yi)
	xi_check = np.array(xi_check)
	yi_check = np.array(yi_check)

	xi2 = []
	yi2 = []
	xi_check2 = []
	yi_check2 = []
	Xvals2 = []
	Yvals2 = []
	xvecs2 = []
	outliersx = []
	outliersy = []

	for i in abs(xi-xi_check):
		if float(i)>1:
			outliersx.append(i)

	for i in abs(yi-yi_check):
		if float(i)>1:
			outliersy.append(i)

	for j in range(len(xi_check)):
		if abs(xi-xi_check)[j] in outliersx or abs(yi-yi_check)[j] in outliersy:
			pass
		else:
			xi2.append(xi[j])
			xi_check2.append(xi_check[j])
			Xvals2.append(Xvals[j])
			yi2.append(yi[j])
			yi_check2.append(yi_check[j])
			Yvals2.append(Yvals[j])
			xvecs2.append(xvecs[j])

	xi = np.array(xi2)
	yi = np.array(yi2)
	xi_check = np.array(xi_check2)
	yi_check = np.array(yi_check2)
	Xvals = np.array(Xvals2)
	Yvals = np.array(Yvals2)
	xvecs = np.array(xvecs2)

	xerror = np.mean(abs(xi-xi_check))
	yerror = np.mean(abs(yi-yi_check))

	print 'Error in x position is {0} pixels'.format(xerror)
	print 'Error in y position is {0} pixels'.format(yerror)

		#convert matched catalogue values to standard coordinates

	Xvecs = []

	for i in range(len(Xvals)):
		Xvec = np.matrix([Xvals[i],Yvals[i],0]).T
		Xvecs.append(Xvec)

	capitalXs = []
	for i in range(len(Xvals)):
		capitalX = np.array((np.linalg.inv(T)*Xvecs[i]).T)
		capitalXs.append(capitalX)

	X_centroid = []
	Y_centroid = []

	for i in range(len(Xvals)):
		X_centroid.append(capitalXs[i][0][0])
		Y_centroid.append(capitalXs[i][0][1])

	X_centroid = np.array(X_centroid)
	Y_centroid = np.array(Y_centroid)

	#print X_centroid

	#construct B matrix
	B = []
	for i in range(len(X_centroid)):
		B.append([m*X_centroid[i],m*Y_centroid[i],1])
	B = np.matrix(B)

	xi = np.matrix(xvecs[:,0])
	yi = np.matrix(xvecs[:,1])

	#calculate x and y constants, as well as reduced chi squared for each
	x_consts = (np.linalg.inv(B.T*B))*B.T*xi
	M = xi - B*x_consts
	x_chi = M.T*M
	dof = len(Xvals)-3
	red_chi_x = x_chi/dof
	#print 'a11 = {0}, \na12 = {1}, \nx0 = {2}, \nwith reduced chi squared = {3}'.format(x_consts[0],x_consts[1],x_consts[2],red_chi_x)
	y_consts = (np.linalg.inv(B.T*B))*B.T*yi
	M = yi - B*y_consts
	y_chi = M.T*M
	dof = len(Yvals)-3
	red_chi_y = y_chi/dof
	#print 'a21 = {0}, \na22 = {1}, \ny0 = {2}, \nwith reduced chi squared = {3}'.format(y_consts[0],y_consts[1],y_consts[2],red_chi_y)

	T = np.matrix([[m*x_consts[0],m*x_consts[1],x_consts[2]],
		[m*y_consts[0],m*y_consts[1],y_consts[2]],
		[0,0,1]])

	#check the transformation by switching catalogue stars with new constants
	xi_check = B*x_consts
	yi_check = B*y_consts

	xi = np.array(xi)
	yi = np.array(yi)
	xi_check = np.array(xi_check)
	yi_check = np.array(yi_check)

	xi2 = []
	yi2 = []
	xi_check2 = []
	yi_check2 = []
	Xvals2 = []
	Yvals2 = []
	xvecs2 = []
	outliersx = []
	outliersy = []

	for i in abs(xi-xi_check):
		if float(i)>0.5:
			outliersx.append(i)

	for i in abs(yi-yi_check):
		if float(i)>0.5:
			outliersy.append(i)

	for j in range(len(xi_check)):
		if abs(xi-xi_check)[j] in outliersx or abs(yi-yi_check)[j] in outliersy:
			pass
		else:
			xi2.append(xi[j])
			xi_check2.append(xi_check[j])
			Xvals2.append(Xvals[j])
			yi2.append(yi[j])
			yi_check2.append(yi_check[j])
			Yvals2.append(Yvals[j])
			xvecs2.append(xvecs[j])

	xi = np.array(xi2)
	yi = np.array(yi2)
	xi_check = np.array(xi_check2)
	yi_check = np.array(yi_check2)
	Xvals = np.array(Xvals2)
	Yvals = np.array(Yvals2)
	xvecs = np.array(xvecs2)

	xerror = np.mean(abs(xi-xi_check))
	errx.append(xerror)
	yerror = np.mean(abs(yi-yi_check))
	erry.append(yerror)

	print 'Error in x position is {0} pixels'.format(xerror)
	print 'Error in y position is {0} pixels'.format(yerror)

	xvals = np.array(x_consts)
	yvals = np.array(y_consts)

	constantsx = np.array([xvals[0][0],xvals[1][0],xvals[2][0]])
	constantsy = np.array([yvals[0][0],yvals[1][0],yvals[2][0]])

	ax.append(constantsx)
	ay.append(constantsy)
