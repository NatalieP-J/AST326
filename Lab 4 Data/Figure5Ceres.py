import numpy as np
import matplotlib.pyplot as plt

#ESTIMATED VALUES
a = 2.947
Omega = 80.65
i = 10.56
e = 0.125
omega = 63.20
tau = 2454833
E = np.arange(0,2*np.pi,0.001)

'''
#TRUE VALUES
a = 2.766
Omega = 80.72
i = 10.61
e = 0.079
omega = 73.12
tau = 2454868
'''

def circular(angles,radius):
	x = radius*np.cos(angles)
	y = radius*np.sin(angles)
	return x,y

def ellipseradius(angles,semimajor,eccentricity):
	return semimajor*(1.0-eccentricity*np.cos(angles))

r = ellipseradius(E,a,e)

plt.plot(circular(E,a)[0],circular(E,a)[1],':')
plt.plot(r*np.cos(E),r*np.sin(E))
plt.axis('equal')
plt.xlim(-3,3)
plt.ylim(-3,3)
plt.show()