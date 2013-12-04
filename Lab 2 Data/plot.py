import matplotlib.pyplot as plt
import numpy as np
import manage as man
from scipy.optimize import leastsq

telescope_fit_night2 = [4.27556985e-01,4.71011978e+02]
lab_fit = [2.35795341e-01,7.43093849e-04,-1.25914456e-06,1.00000000e+00,7.87349037e-10,4.89260494e+02]


#to plot spatial CCD plt.imshow([matrix],origin = 'lower',interpolation='nearest')
#plt.colorbar()

def plot_data(index):
	fnames = ['dark.txt','lamp.txt','neon.txt','overhead.txt','sun.txt','darkneon.txt']
	loc = ['dark','lamp','neon','overhead','sun','neon']
	exper = man.LoadData(fnames[index])
	for m in range(1,len(exper)):
		data = np.genfromtxt('{0}/{1}'.format(loc[index],exper[m][0]), skip_header = 17, skip_footer = 1)
		lam = data[:,0]
		inten = data[:,1]
		print '{0} out of {1}'.format(m,len(exper))
		cont = raw_input('Do you wish to plot (y/n)? ')
		if cont=='y' or cont=='':
			plt.clf()
			plt.plot(lam,inten)
			plt.show()
		if cont=='n':
			break

def average_data(index):
	read_noise = 34.5204291272
	gain = 0.833076727186
	px = np.arange(0,2048)
	p = [366.868875, 0.197670766, -1.529612e-05, -2.024708e-10]
	wv = (guess(p,px))*1e-9
	h = 6.62606957e-34
	c = 3e8
	E = h*c/wv
	fnames = ['dark.txt','lamp.txt','neon.txt','overhead.txt','sun.txt','darkneon.txt']
	loc = ['dark','lamp','neon','overhead','sun','neon']
	exper = man.LoadData(fnames[index])
	dark = man.LoadData(fnames[0])
	intensities = []
	for m in range(1,71):
		data = np.genfromtxt('{0}/{1}'.format(loc[index],exper[m][0]), skip_header = 17, skip_footer = 1)
		darks = np.genfromtxt('{0}/{1}'.format(loc[0],dark[m][0]), skip_header = 17, skip_footer = 1)
		inten = data[:,1]-darks[:,1]-read_noise
		inten/=0.1
		print '{0} out of {1}'.format(m,len(exper))
		intensities.append(inten)
	intensities=np.array(intensities)
	avg = []
	for j in range(0,2048):
		avg.append(np.mean(intensities[:,j]))
	avg = np.array(avg)
	avg /= gain
	for i in range(len(avg)):
		avg[i] *= E[i]
	plt.plot(avg)
	plt.show()
	return avg


def guess(params,independent):
	return params[1]*independent+params[2]*independent**2+params[3]*independent**3+params[0]#*independent**4+params[5]

def process(average,dark,ref,params):
	processed = []
	for i in range(len(average)):
		pro = (average[i]-dark[i])
		processed.append(pro)
	pix = np.arange(0,2048)
	wv = guess(params,pix)
	plt.plot(wv,processed)
	plt.xlim(min(wv),max(wv))
	plt.show()


def maxim_local(data, threshold = 2000):
    maxima = []
    inten = man.IterativeFloatAppend(data,1)
    lam = man.IterativeFloatAppend(data,0)
    for i in range(len(inten)-1):
        if inten[i] > inten[i-1] and inten[i] > inten[i+1] and inten[i] > threshold:
            maxima.append([i,lam[i],inten[i]])
    return maxima

def maxim(data, threshold = 2000):
	inten = data
	maxima = []
	for i in range(len(inten)-1):
		if inten[i] > inten[i-1] and inten[i] > inten[i+1] and inten[i] > threshold:
			maxima.append([i,inten[i]])
	return maxima


'''
hg = average_data(3)
neon = average_data(2)

mh = maxim(hg,threshold=3000)
mn = maxim(neon)

mh = np.array(mh)
mn = np.array(mn)

plt.subplot(211)
plt.plot(hg)
plt.xlim(0,2048)
plt.plot(mh[:,0],mh[:,1],'o')
plt.title('Overhead Light Emission Lines')
plt.ylabel('Signal (ADU)')
plt.subplot(212)
plt.plot(neon)
plt.xlim(0,2048)
plt.title('Neon Emission Lines')
plt.xlabel('Pixel')
plt.ylabel('Signal (ADU)')
plt.plot(mn[:,0],mn[:,1],'o')
plt.show()

'''

'''

neon = []
dark = []
for i in range(70,131,1):
	n = np.loadtxt('AST325_Fall2013/Neon_30sec_01_spectra{0}.dat'.format(i))
	d = np.loadtxt('AST325_Fall2013/Dark_30sec_01_spectra{0}.dat'.format(i))
	dark.append(d)
	neon.append(n-d)

stds = []
avgd=[]
for i in range(len(neon[0])):
	std = []
	sum_avgd = 0
	for j in range(len(neon)):
		sum_avgd += neon[j][i]
		std.append(neon[j][i])
	stds.append(np.std(std)/np.sqrt(len(neon[0])))
	avg = sum_avgd/len(neon)
	avgd.append(avg)

m = maxim(avgd,60)

pixel = []
inten = []
for i in range(len(m)):
	pixel.append(m[i][0])
	inten.append(m[i][1])
	print '{0}'.format(m[i])

plt.plot(pixel,inten,'o')
plt.plot(avgd)
plt.title('Identified maxima of the dark subtracted average neon spectrum')
plt.ylabel('Intensity (ADU)')
plt.xlabel('Pixel')
plt.xlim(0,765)
#plt.show()

'''
neon = man.LoadDataTab('comparisonpixels.dat')
hg = man.LoadData('Hg_pixel_theory.dat')

n = man.IterativeFloatAppend(neon,0)
npixels = man.IterativeIntAppend(neon,1)

h = man.IterativeFloatAppend(hg,0)
hpixels = man.IterativeIntAppend(hg,1)

wavelength = np.array([])
pixels = np.array([])

for i in range(len(h)):
	wavelength = np.append(wavelength,h[i])
	pixels = np.append(pixels,hpixels[i])
for i in range(len(n)):
	wavelength = np.append(wavelength,n[i])
	pixels = np.append(pixels,npixels[i])

smooth_pixel = np.arange(0,2048,1)

def residuals(params,independent,theory):
	return theory - guess(params,independent)

def linguess(params,independent):
	return params[0]*independent+params[1]


def linresiduals(params,independent,theory):
	return theory - linguess(params,independent)

params0 = [1,1,1,1,1,0]

pixels = []
wavelength = []

data = np.loadtxt('telescope_neon.dat')
for i in range(len(data)):
	pixels.append(data[i][0])
	wavelength.append(data[i][1])
pixels = np.array(pixels)
wavelength = np.array(wavelength)
wavelength = wavelength/10
smooth_pixel = np.arange(0,765)

p, cov_x, info, mesg, success = leastsq(residuals,params0,args = (pixels,wavelength),full_output = 1)

f, cov_x, info, mesg, success = leastsq(linresiduals,params0,args = (pixels,wavelength),full_output = 1)

#guess an error in the index of refraction
variance_n=0.1

#calculate goodness of fit parameters
chi_sqr = []
for i in range(len(pixels)):
    resid = (wavelength[i]-guess(p,pixels)[i])**2
    chi = resid/(variance_n**2)
    chi_sqr.append(chi)

chi_squared = sum(chi_sqr)

print 'Chi squared = {0}'.format(chi_squared)

#degrees of freedom
num_data_pts = len(pixels)
num_params = len(p)

dof = num_data_pts - num_params

print 'Data set has {0} degrees of freedom'.format(dof)

#reduced chi squared
red_chi_sqr = chi_squared/dof

print 'Reduced chi squared = {0}'.format(red_chi_sqr)

#residual variance
res_var = np.sqrt(chi_squared/dof)

print 'Residual variance = {0}'.format(res_var)

plt.clf()
plt.subplot(211)
plt.plot(pixels,wavelength,'ob',label='Neon')
plt.plot(smooth_pixel,linguess(f,smooth_pixel),'r',label = 'Linear fit')
plt.plot(smooth_pixel,guess(p,smooth_pixel),'g',label = 'Fit polynomial')
plt.title('Wavelength vs pixel for the telescope spectrometer')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.ylabel('Wavelength (nm)')
plt.legend(loc = 'best')
plt.xlim(0,765)

plt.subplot(212)
plt.plot(pixels,linresiduals(f,pixels,wavelength),'ob',label = 'Mercury residuals')
plt.title('Linear fit residuals')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.ylabel('Wavelength (nm)')
#plt.legend(loc = 'best')
plt.xlim(0,765)
plt.xlabel('Pixel')
plt.subplot(313)

plt.show()
'''

params0 = [1,1,1,1,1,0]

p, cov_x, info, mesg, success = leastsq(residuals,params0,args = (pixels,wavelength),full_output = 1)

uncertainty = 1e-8

print 'cov_x={0}'.format(cov_x)

# Finding chi squared from the residuals and the errors
chi_sqr = sum((pixels - guess(p,pixels))**2 / (uncertainty**2))
print "Chi squared is", "%.4g" % chi_sqr

dof = len(pixels)-len(p)

# Finding the reduced chi squared
red_chi_sqr = chi_sqr/dof
print "Reduced chi squared is", "%.4g" % red_chi_sqr

# And finding the residual variance
res_var = np.sqrt(chi_sqr/dof)
print "Reduced variance is","%.4g" % res_var

# Multiply to get parameter errors
#covs = cov_x*res_var

# And printing the parameter errors (the diagonal elements of the matrix covs)
#for i in range(len(plsq)):
    #print 'The error in p[{0}] is {1}'.format(i,np.sqrt(covs[i,i]))

r = residuals(p,pixels,wavelength)

print residuals(p,pixels,wavelength)

f, cov_x, info, mesg, success = leastsq(linresiduals,params0,args = (pixels,wavelength),full_output = 1)

hpixels = np.array(hpixels)
npixels = np.array(npixels)
n = np.array(n)
h = np.array(h)

plt.clf()
plt.subplot(311)
plt.plot(hpixels,h,'ob',label='Mercury')
plt.plot(npixels,n,'og',label='Neon')
plt.plot(smooth_pixel,linguess(f,smooth_pixel),':',label = 'Linear fit')
plt.plot(smooth_pixel,guess(p,smooth_pixel),'r',label = 'Fit polynomial')
plt.title('Wavelength vs pixel for the lab spectrometer')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.ylabel('Wavelength (nm)')
plt.legend(loc = 'best')
plt.xlim(0,2048)

plt.subplot(312)
plt.plot(hpixels,linresiduals(f,hpixels,h),'ob',label = 'Mercury residuals')
plt.plot(npixels,linresiduals(f,npixels,n),'og',label = 'Neon residuals')
plt.title('Linear fit residuals')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.ylabel('Wavelength (nm)')
plt.legend(loc = 'best')
plt.xlim(0,2048)


plt.subplot(313)
plt.plot(hpixels,residuals(p,hpixels,h),'ob',label = 'Mercury residuals')
plt.plot(npixels,residuals(p,npixels,n),'og',label = 'Neon residuals')
plt.title('Polynomial fit residuals')
plt.xlabel('Pixel')
plt.ylabel('Wavelength (nm)')
plt.legend(loc = 'best')
plt.xlim(0,2048)

plt.tight_layout()
plt.show()

'''