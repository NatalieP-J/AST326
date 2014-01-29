import numpy as np
import matplotlib.pyplot as plt

#ESTIMATED VALUES
a = 2.947
Omega = 80.65
i = 10.56
e = 0.125
omega = 63.20
tau = 2454833


#TRUE VALUES
a = 2.766
Omega = 80.72
i = 10.61
e = 0.079
omega = 73.12
tau = 2454868


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
k = 0.01720209895
n = np.sqrt(k**2/a**3)
period = (2*np.pi)/n

#Generate time intervals
start = 2454702.5
JulianDay = np.arange(start,start+period+500,1)
J = JulianDay - 2450000

M = meananomalytime(n,tau,start)
E0 = M

E = []
E.append(E0)

for i in range(1,len(JulianDay)):
	delE = (meananomalytime(n,tau,JulianDay[i]) - meananomalyE(E[i-1],e))/(1-e*np.cos(E[i-1]))
	En = E[i-1] + delE
	E.append(En)

E = np.array(E)
M = meananomalytime(n,tau,JulianDay)

v = trueanomaly(e,E)
r = ellipseradius(E,a,e)

theta = v+(omega*(np.pi/180))

perihelion = []
vcheck = []
for i in range(len(v)):
	if np.round(v[i],2) == 0:
		vcheck.append(v[i])

for i in range(len(v)):
	if v[i] == min(vcheck):
		perihelion.append(r[i]*np.cos(theta[i]))
		perihelion.append(r[i]*np.sin(theta[i]))

#FIGURE 3

plt.figure()
plt.plot(J,E,'k',label='Eccentric Anomaly')
plt.plot(J,M,'k:',label='Mean Anomaly')
plt.ylim(0,8)
plt.xlim(4500,7000)
plt.xlabel('Julian Day - 2450000')
plt.ylabel('Anomaly [radians]')
plt.title('Anomaly Over Time')
plt.legend(loc='best')
plt.show()

#FIGURE 4

plt.figure()
plt.subplot(211)
plt.plot(J,r,'k')
plt.ylim(2.4,3.0)
plt.xlim(4500,7000)
plt.xlabel('Julian Day - 2450000')
plt.ylabel('r [AU]')
plt.title('Orbital Separation Over Time')
plt.subplot(212)
plt.plot(J,v,'k')
plt.ylim(-4,4)
plt.xlim(4500,7000)
plt.xlabel('Julian Day - 2450000')
plt.ylabel('v [radians]')
plt.title('True Anomaly Over Time')
plt.tight_layout()
plt.show()

#FIGURE 5

plt.figure()
plt.plot(circular(theta,a)[0],circular(theta,a)[1],':',label = 'Circular orbit radius = {0}'.format(a))
plt.plot(r*np.cos(theta),r*np.sin(theta),'k',label = 'Calculated orbit')
plt.plot(0,0,'k+')
plt.plot(perihelion[0],perihelion[1],'ko')
plt.axis('equal')
plt.xlim(-3,3)
plt.ylim(-3,3)
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.legend(loc = 'best')
plt.title('Position of Ceres in its Orbital Plane')
plt.show()
