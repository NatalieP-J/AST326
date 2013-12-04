import numpy as np
import matplotlib.pyplot as plt
from plot import *
from scipy.optimize import leastsq

k = 1.3806488e-23
h = 6.62606957e-34
c = 3e8
b = 2.8977685e-3

def pix2wv(pix,params):
	return (params[0]*pix+params[1])*1e-9

def blackbody(params,x):
	m = 2*h*(c**2)/x**5
	d = np.exp((h*c)/(x*k*params[0]))-1
	b = m*(1./d)-params[1]
	return b

def residuals(temp,x,y):
	return y-blackbody(temp,x)

def absorption(data,threshold,wv):
	minima=[]
	mod_min=[]
	for i in range(1,len(data)-1):
		if data[i]<data[i-1] and data[i]<data[i+1]:
			minima.append([wv[i],data[i]])
	minima=np.array(minima)
	for k in range(len(data)):
		if data[k] in minima[:,1]:
			avg = np.mean(data[k-10:k-1])
			avg += np.mean(data[k+1:k+10])
			avg = avg/2
			if abs(data[k]-avg) >= threshold:
				mod_min.append([wv[k],data[k]])
	return np.array(mod_min)

def emission(data,threshold,wv):
	maxima=[]
	mod_max=[]
	for i in range(1,len(data)-1):
		if data[i]>data[i-1] and data[i]>data[i+1]:
			maxima.append([wv[i],data[i]])
	maxima=np.array(maxima)
	for k in range(len(data)):
		if data[k] in maxima[:,1]:
			avg = np.mean(data[k-10:k-1])
			avg += np.mean(data[k+1:k+10])
			avg = avg/2
			if abs(data[k]-avg) >= threshold:
				mod_max.append([wv[k],data[k]])
	return np.array(mod_max)

neon = []
neon_dark =[]
for i in range(24,124):
	neon.append(np.loadtxt('Neon1_2min_spectra{0}.dat'.format(i)))
	neon_dark.append(np.loadtxt('Neon1_dark_2min_spectra{0}.dat'.format(i)))
neon = np.array(neon)
neon_dark = np.array(neon_dark)

neon -= neon_dark

neon_avg = []
for i in range(0,765):
	neon_avg.append(np.mean(neon[:,i]))

neon_avg = np.array(neon_avg)

#with a leastsq fit:
p = [0.42743091,383.59310499]

pixels = np.arange(0,765)
wv = pix2wv(pixels,p)

vega = []
vflat = []
vega_dark = []
vflat_dark = []
for i in range(209,310):
	vega.append(np.loadtxt('vega03_1min_spectra{0}.dat'.format(i)))
	vega_dark.append(np.loadtxt('vega03_dark_spectra{0}.dat'.format(i)))
	vflat.append(np.loadtxt('Flat1_300sec_spectra{0}.dat'.format(i)))
	vflat_dark.append(np.loadtxt('Flat1_dark_300sec_spectra{0}.dat'.format(i)))

vega = np.array(vega)
vega = vega/30
vflat = np.array(vflat)
vflat = vflat/300
vega_dark = np.array(vega_dark)
vflat_dark = np.array(vflat_dark)

vega -= vega_dark/30
vflat -= vflat_dark/300


vega_avg = []
vflat_avg = []

for i in range(0,765):
	vega_avg.append(np.mean(vega[:,i]))
	vflat_avg.append(np.mean(vflat[:,i]))

vega_avg = np.array(vega_avg)
vflat_avg = np.array(vflat_avg)

#temp of lamp
m = max(vflat_avg)
for i in range(len(vflat_avg)):
	if vflat_avg[i]==m:
		max_wv = wv[i]

temp = b/max_wv
planck = blackbody([temp,0],wv)

processed = (vega_avg/vflat_avg)*planck

vm = max(processed)
for i in range(len(processed)):
	if processed[i]==vm:
		vmax_wv = wv[i]

vtemp = b/vmax_wv

print 'Vega temperature is {0}K'.format(vtemp)

params = [vtemp,2.5e13]

p = leastsq(residuals,params,args=(wv,processed))

minima = absorption(processed,3.5e12,wv)

'''
plt.figure()
plt.subplot(211)
plt.plot(vega[54])
plt.subplot(212)
plt.plot(vega_avg)
#plt.show()

'''

plt.figure()
plt.subplot(211)
plt.plot(wv,vega_avg)
plt.xlim(min(wv),max(wv))
#plt.xlabel('Wavelength (nm)')
plt.ylabel('Signal (ADU)')
plt.title('Average Dark Subracted Vega Spectrum')
plt.subplot(212)
plt.xlim(min(wv),max(wv))
plt.xlabel('Wavelength (m)')
plt.ylabel(r'Intensity $(W\cdot m^-3)$')
plt.title('Processed Vega Spectrum')
plt.plot(wv,processed)
plt.plot(minima[:,0],minima[:,1],'o')
plt.text(6.65e-7,1.7e13,r'H$\alpha$')
plt.text(4.95e-7,2.5e13,r'H$\beta$')
plt.text(4.45e-7,1.5e13,r'H$\gamma$')
plt.text(4.20e-7,1e13,r'H$\delta$')
plt.show()


'''
erif1=[]
e1f=[]
for i in range(270,300):
	erif1.append(np.loadtxt('Scheat_2min_01_spectra{0}.dat'.format(i)))
	e1f.append(np.loadtxt('Flat1_300sec_spectra{0}.dat'.format(i)))

erif1=np.array(erif1)
erif1=erif1/120
e1f=np.array(e1f)
e1f=e1f/300

erif1_avg = []
e1f_avg = []

for i in range(765):
	erif1_avg.append(np.mean(erif1[:,i]))
	e1f_avg.append(np.mean(e1f[:,i]))

erif1_avg = np.array(erif1_avg)
e1f_avg = np.array(e1f_avg)

p = (erif1_avg/e1f_avg)*blackbody([vtemp,0],wv)
maxima = emission(p,3e11,wv)
minima = absorption(p,3e11,wv)


p_spec = []
for i in range(len(p)):
	if p[i] not in maxima:
		p_spec.append(p[i])

p_spec = np.array(p_spec)

vm = max(p_spec)
for i in range(len(p_spec)):
	if p_spec[i]==vm:
		vmax_wv = wv[i]

vtemp = b/vmax_wv

print 'Scheat temperature is {0}K'.format(vtemp)

plt.subplot(411)
plt.plot(wv,p)
plt.xlim(min(wv),max(wv))
if len(maxima) != 0:
	plt.plot(maxima[:,0],maxima[:,1],'o')
if len(minima) != 0:
	plt.plot(minima[:,0],minima[:,1],'o')
#plt.xlabel('Wavelength (m)')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.ylabel('Signal (ADU)')
plt.title('Scheat')

'''
erif2 = []
e2f = []
for i in range(250,269):
	erif2.append(np.loadtxt('Erif02_spectra{0}.dat'.format(i)))
	e2f.append(np.loadtxt('Flat1_300sec_spectra{0}.dat'.format(i)))

erif2 = np.array(erif2)
erif2 = erif2/60
e2f = np.array(e2f)
e2f = e2f/300

erif2_avg = []
e2f_avg = []

for i in range(765):
	erif2_avg.append(np.mean(erif2[:,i]))
	e2f_avg.append(np.mean(e2f[:,i]))

erif2_avg = np.array(erif2_avg)
e2f_avg = np.array(e2f_avg)

p = (erif2_avg/e2f_avg)*blackbody([vtemp,0],wv)
maxima = emission(p,3e12,wv)
minima = absorption(p,2.1e12,wv)

p_spec = []
for i in range(len(p)):
	if p[i] not in maxima:
		p_spec.append(p[i])

p_spec = np.array(p_spec)

vm = max(p_spec)
for i in range(len(p_spec)):
	if p_spec[i]==vm:
		vmax_wv = wv[i]

vtemp = b/vmax_wv

print 'Enif temperature is {0}K'.format(vtemp)

plt.subplot(311)
plt.plot(wv,p)
plt.xlim(min(wv),max(wv))
if len(maxima) != 0:
	plt.plot(maxima[:,0],maxima[:,1],'o')
if len(minima) != 0:
	plt.plot(minima[:,0],minima[:,1],'o')
#plt.xlabel('Wavelength (m)')
plt.ylabel(r'Intensity $(W\cdot m^-3)$')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.title('Enif 2')

erif2 = []
e2f = []
for i in range(290,311):
	erif2.append(np.loadtxt('albero_01_spectra{0}.dat'.format(i)))
	e2f.append(np.loadtxt('Flat1_300sec_spectra{0}.dat'.format(i)))

erif2 = np.array(erif2)
erif2 = erif2/60
e2f = np.array(e2f)
e2f = e2f/300

erif2_avg = []
e2f_avg = []

for i in range(765):
	erif2_avg.append(np.mean(erif2[:,i]))
	e2f_avg.append(np.mean(e2f[:,i]))

erif2_avg = np.array(erif2_avg)
e2f_avg = np.array(e2f_avg)

p = (erif2_avg/e2f_avg)*blackbody([vtemp,0],wv)
maxima = emission(p,1e11,wv)
minima = absorption(p,1e11,wv)

p_spec = []
for i in range(len(p)):
	if p[i] not in maxima:
		p_spec.append(p[i])

p_spec = np.array(p_spec)

vm = max(p_spec)
for i in range(len(p_spec)):
	if p_spec[i]==vm:
		vmax_wv = wv[i]
print vmax_wv
vtemp = b/vmax_wv

print 'Albireo 01 temperature is {0}K'.format(vtemp)

plt.subplot(312)
plt.plot(wv,p)
plt.xlim(min(wv),max(wv))
if len(maxima) != 0:
	plt.plot(maxima[:,0],maxima[:,1],'o')
if len(minima) != 0:
	plt.plot(minima[:,0],minima[:,1],'o')
#plt.xlabel('Wavelength (m)')
plt.ylabel(r'Intensity $(W\cdot m^-3)$')
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off') 
plt.title('Albireo 1')

erif2 = []
e2f = []
for i in range(245,270):
	erif2.append(np.loadtxt('albero_02_spectra{0}.dat'.format(i)))
	e2f.append(np.loadtxt('Flat1_300sec_spectra{0}.dat'.format(i)))

erif2 = np.array(erif2)
erif2 = erif2/60
e2f = np.array(e2f)
e2f = e2f/300

erif2_avg = []
e2f_avg = []

for i in range(765):
	erif2_avg.append(np.mean(erif2[:,i]))
	e2f_avg.append(np.mean(e2f[:,i]))

erif2_avg = np.array(erif2_avg)
e2f_avg = np.array(e2f_avg)

p = (erif2_avg/e2f_avg)*blackbody([vtemp,0],wv)
maxima = emission(p,1e12,wv)
minima = absorption(p,1e12,wv)

p_spec = []
for i in range(len(p)):
	if p[i] not in maxima:
		p_spec.append(p[i])

p_spec = np.array(p_spec)

vm = max(p_spec)
for i in range(len(p_spec)):
	if p_spec[i]==vm:
		vmax_wv = wv[i]

vtemp = b/vmax_wv
print vmax_wv
print 'Albireo 02 temperature is {0}K'.format(vtemp)

plt.subplot(313)
plt.plot(wv,p)
plt.xlim(min(wv),max(wv))
if len(maxima) != 0:
	plt.plot(maxima[:,0],maxima[:,1],'o')
if len(minima) != 0:
	plt.plot(minima[:,0],minima[:,1],'o')
plt.xlabel('Wavelength (m)')
plt.ylabel(r'Intensity $(W\cdot m^-3)$')
plt.title('Albireo 2')

#plt.tight_layout()
plt.show()



