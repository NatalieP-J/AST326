import numpy as np
import matplotlib.pyplot as plt

k = 1.3806488e-23
h = 6.62606957e-34
c = 3e8

telescope_fit = [4.27556985e-01,4.71011978e+02]

def pix2wv(range,params):
	pix = np.arange(0,range,1)
	wv = (params[0]*pix+params[1])*1e-9
	return np.array(wv)

def blackbody(temp,x):
	m = 2*h*(c**2)/x**5
	d = np.exp((h*c)/(x*k*temp))-1
	b = m*(1./d)
	return b

wv = pix2wv(765,telescope_fit)

name = 'AST325_Fall2013'

vega = []
vdark30 = []
vdark2 = []
vflat = []
scheat = []
sdark2 = []
sflat = []

for i in range(103,204):
	vega.append(np.loadtxt('{0}/Vega_30sec_01_spectra{1}.dat'.format(name,i)))
	vdark30.append(np.loadtxt('{0}/Dark_30sec_01_spectra{1}.dat'.format(name,i)))
	vdark2.append(np.loadtxt('{0}/Dark_2min_02_spectra{1}.dat'.format(name,i)))
	vflat.append(np.loadtxt('{0}/Flat_2m_01_spectra{1}.dat'.format(name,i)))

for i in range(124,325):
	scheat.append(np.loadtxt('{0}/Scheat_2min_01_spectra{1}.dat'.format(name,i)))
	sdark2.append(np.loadtxt('{0}/Dark_2min_02_spectra{1}.dat'.format(name,i)))
	sflat.append(np.loadtxt('{0}/Flat_2m_01_spectra{1}.dat'.format(name,i)))

vega = np.array(vega)
vdark30 = np.array(vdark30)
vdark2 = np.array(vdark2)
vflat = np.array(vflat)
scheat = np.array(scheat)
sdark2 = np.array(sdark2)
sflat = np.array(sflat)

vega -= vdark30
vflat -= vdark2
scheat -= sdark2
sflat -= sdark2

vegaAvg = []
vflatAvg = []
scheatAvg = []
sflatAvg = []
for i in range(0,765):
	vegaAvg.append(np.mean(vega[:,i]))
	vflatAvg.append(np.mean(vflat[:,i]))
	scheatAvg.append(np.mean(scheat[:,i]))
	sflatAvg.append(np.mean(sflat[:,i]))
vegaAvg = np.array(vegaAvg)
vflatAvg = np.array(vflatAvg)
scheatAvg = np.array(scheatAvg)
sflatAvg = np.array(sflatAvg)




