import numpy as np
import matplotlib.pyplot as plt
from plot import *

neon = []
neon_dark =[]
for i in range(24,124):
	neon.append(np.loadtxt('Neon1_2min_spectra{0}.dat'.format(i)))
	neon_dark.append(np.loadtxt('Neon1_dark_2min_spectra{0}.dat'.format(i)))
neon = np.array(neon)
neon_dark = np.array(neon_dark)
