import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
import urllib as url
import string as str
from scipy.optimize import leastsq

x_consts = np.array([1.00050002e+00,-6.59823982e-03,1.01964337e+03])
y_consts = np.array([6.26632180e-03,9.99831623e-01,1.02834801e+03])
asteroid_x = [1089.5,1090.5,1065,1088.5,1071.5]
asteroid_y = [977.5,1039.5,1008,1039,1025]
actual_ax = [1087.1291575346165,1089.0752263384331,1063.9500813357831,1087.0836146043609,1070.5253099208114]
actual_ay = [992.70536623976943,1008.3996993378165,1040.0641397572322,1008.8906946228163,1022.642015275929]
radian = np.pi/180
errors = [0.755029807194,0.756458254734]
f = 3420
p = 0.018
m = f/p
T = []
t = [m*x_consts[0],m*x_consts[1],x_consts[2]]
T.append(t)
t = [m*y_consts[0],m*y_consts[1],y_consts[2]]
T.append(t)
t = [0,0,1]
T.append(t)
T = np.matrix(T)
Tinv = np.linalg.inv(T)

skies = []
darks = []
flats = []
corrected = []

days = [19,20,22,23,28]

def pix2rad(p0,X,Y):
	ra0,dec0 = p0
	ra = np.arctan(-X/(np.cos(dec0)-Y*np.sin(dec0)))+ra0
	dec = np.arcsin((np.sin(dec0)+Y*np.cos(dec0))/((1+(X**2)+(Y**2))**0.5))
	return ra,dec

righta = []
decl = []
Xasteroid = []
Yasteroid = []
actual_asteroidx = []
actual_asteroidy = []
xerrors = []
yerrors = []

for urania in range(5):
	print 'Working on {0}'.format(urania)
	skyname = '30Urania-S001-R001-C001-r-{0}.fts'.format(urania)
	darkname = 'Dark-S001-R001-C001-B2-{0}.fts'.format(urania)
	flatname = 'AutoFlat-Dusk-r-Bin2-001-{0}.fts'.format(urania)
	s = pf.open(skyname)
	sky_time = s[0].header['exptime']
	sky_date = s[0].header['date']
	sky_hour = s[0].header['time-obs']
	print 'Raw exposure {0}s'.format(sky_time)
	print 'Raw date {0}'.format(sky_date)
	print 'Raw time {0}'.format(sky_hour)
	sky = pf.getdata(skyname)
	sky /= sky_time
	s = pf.open(darkname)
	dark_time = s[0].header['exptime']
	print 'Dark exposure {0}s'.format(dark_time)
	dark = pf.getdata(darkname)
	dark /= dark_time
	s = pf.open(flatname)
	flat_time = s[0].header['exptime']
	print 'Flat exposure {0}s'.format(flat_time)
	flat = pf.getdata(flatname)
	flat /= flat_time
	skies.append(sky)
	darks.append(dark)
	flats.append(flat)
	f = flat-dark
	f = f/np.median(f)
	correct = ((sky-dark)/f)[::-1]
	corrected.append(correct)
	x = np.loadtxt('Jan{0}_Xcentroid.txt'.format(days[urania]))
	y = np.loadtxt('Jan{0}_Ycentroid.txt'.format(days[urania]))
	X = Tinv[0,0]*x+Tinv[0,1]*y+Tinv[0,2]
	Y = Tinv[1,0]*x+Tinv[1,1]*y+Tinv[1,2]
	aX = Tinv[0,0]*asteroid_x[urania]+Tinv[0,1]*asteroid_y[urania]+Tinv[0,2]
	aY = Tinv[1,0]*asteroid_x[urania]+Tinv[1,1]*asteroid_y[urania]+Tinv[1,2]
	actual_aX = Tinv[0,0]*actual_ax[urania]+Tinv[0,1]*actual_ay[urania]+Tinv[0,2]
	actual_aY = Tinv[1,0]*actual_ax[urania]+Tinv[1,1]*actual_ay[urania]+Tinv[1,2]
	xerror = Tinv[0,0]*errors[0]+Tinv[0,1]*errors[1]+Tinv[0,2]
	yerror = Tinv[1,0]*errors[0]+Tinv[1,1]*errors[1]+Tinv[1,2]
	s1 = pf.open(skyname)
	ras = s1[0].header['ra']
	des = s1[0].header['dec']
	radeg = 15*(float(ras[0:2])+float(ras[3:5])/60. + float(ras[6:])/3600.)
	dsgn = np.sign(float(des[0:3]))
	dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.
	ra0 = radeg*radian
	dec0 = dedeg*radian
	ra,dec = pix2rad([ra0,dec0],X,Y)
	Xa,Ya = pix2rad([ra0,dec0],aX,aY)
	actualx,actualy = pix2rad([ra0,dec0],actual_aX,actual_aY)
	checkx,checky = pix2rad([ra0,dec0],(actual_aX+xerror),(actual_aY+yerror))
	Xerror = abs(actualx-checkx)
	Yerror = abs(actualy-checky)
	Xerror /= radian
	Yerror /= radian
	ra /= radian
	dec /= radian
	Xa /= radian
	Ya /= radian
	actualx /= radian
	actualy /= radian
	xerrors.append(Xerror)
	yerrors.append(Yerror)
	actual_asteroidx.append(actualx)
	actual_asteroidy.append(actualy)
	Xasteroid.append(Xa)
	Yasteroid.append(Ya)
	righta.append(ra)
	decl.append(dec)

Xasteroid = np.array(Xasteroid)
Yasteroid = np.array(Yasteroid)
actual_asteroidx = np.array(actual_asteroidx)
actual_asteroidy = np.array(actual_asteroidy)

#GET PROPER MOTION

#time between observations in seconds
time = np.array([87117,176593,81788])
timehr = time/3600.

vx = []
vy = []
for i in range(3):
	speedx = (actual_asteroidx[i+1]-actual_asteroidx[i])/timehr[i]
	vx.append(speedx)
	speedy = (actual_asteroidy[i+1]-actual_asteroidy[i])/timehr[i]
	vy.append(speedy)

actual_asteroidx = np.array(actual_asteroidx)
actual_asteroidy = np.array(actual_asteroidy)
time = np.array([0,87117,(176593+87117),(81788+176593+87117),(421230+81788+176593+87117)])
timehr = time/3600.

def line(p,x):
	m,b = p
	return (m*x)+b

def residuals(p,x,y):
	return y-line(p,x)

p0 = [1e-6,44]

popt, cov_x, info, mesg, success = leastsq(residuals, p0, args = (timehr,actual_asteroidx), full_output = True,maxfev=2000)

#Goodness of fit analysis

#guess an error in the Intensity based on the data.
std_n=0.870389808/3600

#calculate goodness of fit parameters
chi_sqr = []
for i in range(len(actual_asteroidx)):
    resid = (residuals(popt,timehr,actual_asteroidx)[i])**2
    chi = resid/(std_n**2)
    chi_sqr.append(chi)

chi_squared = sum(chi_sqr)

print 'Chi squared = {0}'.format(chi_squared)

#degrees of freedom
num_data_pts = len(actual_asteroidx)
num_params = len(popt)

dof = num_data_pts - num_params

print 'Data set has {0} degrees of freedom'.format(dof)

#reduced chi squared
red_chi_sqr = chi_squared/dof

print 'Reduced chi squared = {0}'.format(red_chi_sqr)

#residual variance
res_var = np.sqrt(chi_squared/dof)

print 'Residual variance = {0}'.format(res_var)

# Multiply to get parameter errors
covs = cov_x*res_var

# And printing the parameter errors (the diagonal elements of the matrix covs)
p_round = []
for i in range(len(popt)):
    print 'The error in p[{0}] is {1}'.format(i,np.sqrt(covs[i,i]))

t = np.delete(timehr,len(timehr)-1)
actx = np.delete(actual_asteroidx,len(actual_asteroidx)-1)
acty = np.delete(actual_asteroidy,len(actual_asteroidy)-1)

ropt, cov_x, info, mesg, success = leastsq(residuals, p0, args = (t,actx), full_output = True,maxfev=2000)

p0 = [1e-7,19]

kopt, cov_x, info, mesg, success = leastsq(residuals, p0, args = (timehr,actual_asteroidy), full_output = True,maxfev=2000)

mopt, cov_x, info, mesg, success = leastsq(residuals, p0, args = (t,acty), full_output = True,maxfev=2000)

smooth_timehr = np.arange(-44,timehr[4]+100,0.1)

plt.figure()
plt.subplot(211)
#plt.plot(timehr,actual_asteroidx,'o')
plt.errorbar(timehr,actual_asteroidx,yerr = 0.890389808/3600,fmt='o')
plt.plot(smooth_timehr,line(popt,smooth_timehr),'m',label='ra proper motion: 0.65 \'/hr')
plt.plot(smooth_timehr,line(ropt,smooth_timehr),'g',label='ra proper motion: 0.60 \'/hr')
plt.xlim(min(timehr),max(timehr))
plt.ylabel(r'Right Ascension($^o$)',fontsize=10)
plt.xlabel('Time (hour)',fontsize=10)
plt.title('Motion in Right Ascension',fontsize=11)
plt.legend(loc='best')
plt.subplot(212)
#plt.plot(timehr,actual_asteroidy,'o')
plt.errorbar(timehr,actual_asteroidy,yerr = 0.8173101414/3600,fmt='o')
plt.plot(smooth_timehr,line(kopt,smooth_timehr),'m',label='dec proper motion: 0.11 \'/hr')
plt.plot(smooth_timehr,line(mopt,smooth_timehr),'g',label='dec proper motion: 0.10 \'/hr')
plt.xlim(min(timehr),max(timehr))
plt.ylabel(r'Declination($^o$)',fontsize=10)
plt.xlabel('Time (hour)',fontsize=10)
plt.title('Motion in Declination',fontsize=11)
plt.legend(loc='best')
plt.suptitle('Proper Motion of 30 Urania')
plt.show()


cs = ['b','g','c','m','y']
plt.figure()
for i in range(len(righta)):
	plt.errorbar(righta[i],decl[i],label = 'Jan {0}'.format(days[i]),yerr=0.8173101414/3600,xerr = 0.890389808/3600,fmt='{0}.'.format(cs[i]))
plt.errorbar(actual_asteroidx,actual_asteroidy,fmt='rd',yerr=0.8173101414/3600,xerr = 0.890389808/3600)
plt.plot(line(popt,smooth_timehr),line(kopt,smooth_timehr),'r',label='Fitted asteroid motion')
#plt.plot(line(ropt,smooth_timehr),line(mopt,smooth_timehr),'g',label='Fitted asteroid motion')
plt.xlim(min(righta[0]),max(righta[4]))
plt.ylim(min(decl[0]),max(decl[4]))
plt.xlabel(r'Right Ascension ($^o$)')
plt.ylabel(r'Declination ($^o$)')
plt.title('Path of 30 Urania Across the Sky')
plt.legend(loc='best')
plt.show()

