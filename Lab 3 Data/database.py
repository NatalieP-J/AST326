import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
import urllib as url
import string as str
from scipy.optimize import leastsq

#OUTPUT
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

galaxy = 'NGC7331-S001-R001-C001-r.fts' #exposure time 240s
dark = 'Dark-S001-R003-C003-B2.fts'
flat = 'combflatr.fits'
d = pf.getdata(dark)
dhdr = pf.getheader(dark)
x = pf.getdata(galaxy)
f = pf.getdata(flat)
fnew = f/np.median(f)
correct = ((x-d)/fnew)[::-1]

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

galaxy = 'NGC7331-S001-R001-C001-r.fts'
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
x = np.loadtxt('Xcentroid.txt')
y = np.loadtxt('Ycentroid.txt')

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
		if R < 15:
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
plt.show()

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
plt.show()

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
print 'a11 = {0}, \na12 = {1}, \nx0 = {2}, \nwith reduced chi squared = {3}'.format(x_consts[0],x_consts[1],x_consts[2],red_chi_x)
y_consts = (np.linalg.inv(B.T*B))*B.T*yi
M = yi - B*y_consts
y_chi = M.T*M
dof = len(Yvals)-3
red_chi_y = y_chi/dof
print 'a21 = {0}, \na22 = {1}, \ny0 = {2}, \nwith reduced chi squared = {3}'.format(y_consts[0],y_consts[1],y_consts[2],red_chi_y)

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
plt.show()

xerror = np.mean(abs(xi-xi_check))
yerror = np.mean(abs(yi-yi_check))

print 'Error in x position is {0} pixels'.format(xerror)
print 'Error in y position is {0} pixels'.format(yerror)

plt.figure()
plt.plot(x,y,'b.',label = 'Centroids')
plt.plot(xi_check,yi_check,'.r',label='Transformed Catalogue Stars')
plt.xlabel('x [pixel]')
plt.ylabel('y [pixel]')
plt.title('Overplot of Centroids and Galaxy Stars With New Transformation Constants')
plt.ylim(0,2048)
plt.xlim(0,2048)
plt.legend(loc='best')
plt.show()
