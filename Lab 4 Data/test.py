import numpy as np
import matplotlib.pyplot as plt

#ESTIMATED VALUES
a = 2.947
Omega = 80.65
i = 10.56
e = 0.125
omega = 63.20
tau = 2454833

#Constants
Msun = 2e30
G = 6.67e-11
k = np.sqrt(G*Msun)
k = 0.01720209895
n = np.sqrt(k**2/a**3)
period = k*2*np.pi*a**1.5


#Generate time intervals
start = 2454702.5
JulianDay = np.arange(start,start+period+500,1)