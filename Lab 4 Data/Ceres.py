import numpy as np
import matplotlib.pyplot as plt

#ESTIMATED VALUES
a = 2.947
Omega = 80.65
i = 10.56
e = 0.125
omega = 63.20
tau = 2454833

'''
#TRUE VALUES
a = 2.766
Omega = 80.72
i = 10.61
e = 0.079
omega = 73.12
tau = 2454868
'''



def ellipseradius(angles,semimajor,eccentricity):
	return semimajor*(1.0-eccentricity*np.cos(angles))

def trueanomaly(eccentricity,angles):
	return 2*np.arctan(np.sqrt((1+eccentricity)/(1-eccentricity))*np.tan(angles/2))

def meananomalytime(n,tau,t):
	return n*(t-tau)

def meananomalyE(E,e):
	return E-e*np.sin(E)

def circular(angles,radius):
	x = radius*np.cos(angles)
	y = radius*np.sin(angles)
	return x,y


#Constants
Msun = 2e30
G = 6.67e-11
k = np.sqrt(G*Msun)
n = np.sqrt(k**2/a**3)
period = k*2*np.pi*a**1.5

#Generate time intervals
start = 2454702.5
JulianDay = np.arange(start,start+period+500,1)

M = meananomalytime(n,tau,start)
print 'mean anomaly time'
E0 = M

E = []
E.append(E0)

for i in range(1,len(JulianDay)):
	delE = (meananomalytime(n,tau,JulianDay[i]) - meananomalyE(E[i-1],e))/(1-e*np.cos(E[i-1]))
	En = E[i-1] + delE
	E.append(En)

print 'eccentric anomaly'
E = np.array(E)
M = meananomalytime(n,tau,JulianDay)
print 'mean anomaly'

v = trueanomaly(e,E)
print 'true anomaly'
r = ellipseradius(E,a,e)
print 'radius'

theta = v+(omega*(np.pi/180))

'''
#FIGURE 3

plt.figure()
plt.plot(JulianDay,E,label='Eccentric Anomaly')
plt.plot(JulianDay,M,label='Mean Anomaly')
plt.show()

#FIGURE 4

plt.figure()
plt.subplot(211)
plt.plot(Julian,r)
plt.subplot(212)
plt.plot(Julian,v)
plt.show()

#FIGURE 5

plt.figure()
plt.plot(circular(theta,a)[0],circular(theta,a)[1],':')
plt.plot(r*np.cos(theta),r*np.sin(theta))
plt.axis('equal')
plt.xlim(-3,3)
plt.ylim(-3,3)
plt.show()
'''