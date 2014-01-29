import numpy as np
import matplotlib.pyplot as plt


#DEFINITIONS AND CONSTANTS
cos = np.cos
sin = np.sin
rad = np.pi/180

k = 0.01720209895

#FUNCTIONS 

def eclipticToUnit(lam,beta):
	x = cos(lam)*cos(beta)
	y = sin(lam)*cos(beta)
	z = sin(beta)
	return x,y,z

def rho(k,R,s,sdot,sdouble,r):
	numerator = np.dot(sdot,np.cross(R,s))
	denominator = np.dot(sdot,np.cross(sdouble,s))
	Rmag = np.linalg.norm(R)
	constants = (k**2)*((1./Rmag**3)-(1./r**3))
	return constants*(numerator/denominator)

def rhodot(k,R,s,sdot,sdouble,r):
	numerator = np.dot(sdouble,np.cross(R,s))
	denominator = np.dot(sdouble,np.cross(sdot,s))
	Rmag = np.linalg.norm(R)
	constants = ((k**2)/2)*((1./Rmag**3)-(1./r**3))
	return constants*(numerator/denominator)

def radius(rhoval,R,s):
	Rmag = np.linalg.norm(R)
	return np.sqrt((rhoval**2)+(Rmag**2)+2*rhoval*np.dot(R,s))

def trueanomaly(eccentricity,angles):
	return 2*np.arctan(np.sqrt((1+eccentricity)/(1-eccentricity))*np.tan(angles/2))

def meananomalyE(E,e):
	return E-e*np.sin(E)

#MEASUREMENT DATES

Julian = np.array([2454702.5,2454703.5,2454704.5])

#CERES POSITION (ecliptic - radians)
ceres = np.array([[121.7592648*rad,4.0625653*rad],[122.1865441*rad,4.0992581*rad],[122.6133849*rad,4.1361592*rad]])
lam = ceres[:,0]
beta = ceres[:,1]

#EARTH POSITION (ecliptic - AU)
earth = np.array([[0.8849686471,-0.4888489729,4.466373306e-06],[0.8928865393,-0.4737871683,4.402701086e-6],[0.9005490495,-0.4585878955,4.48380158e-6]])
X = earth[:,0]
Y = earth[:,1]
Z = earth[:,2]

#FINDING UNIT VECTOR s

x,y,z = eclipticToUnit(lam,beta)
s1 = np.array([x[0],y[0],z[0]])
s2 = np.array([x[1],y[1],z[1]])
s3 = np.array([x[2],y[2],z[2]])

#GENERATING DERIVATIVES OF UNIT VECTOR s

tau1 = Julian[1]-Julian[0]
tau3 = Julian[2]-Julian[1]

s2dot = ((tau3/(tau1*(tau1+tau3)))*(s2-s1)) + ((tau1/(tau3*(tau1+tau3)))*(s3-s2))
s2double = ((2/(tau3*(tau1+tau3)))*(s3-s2)) - ((2/(tau1*(tau1+tau3)))*(s2-s1))

#r,rho AND rhodot

r0 = 2.5 #AU
radiuslist = []
rholist = []
radiuslist.append(r0)

R2 = earth[1]

for i in range(10):
	p = rho(k,R2,s2,s2dot,s2double,radiuslist[i])
	rholist.append(p)
	rval = radius(rholist[i],R2,s2)
	radiuslist.append(rval)

rhovel = rhodot(k,R2,s2,s2dot,s2double,radiuslist[len(rholist)-1])

#print 'r = {0}, rho = {1}, rhodot = {2}'.format(radiuslist[9],rholist[9],rhovel)

plt.figure()
plt.subplot(211)
plt.plot(rholist,'.')
plt.title('distance from Earth iterative solution')
plt.ylabel(r'$\rho$')
plt.xlabel('Iteration')
plt.subplot(212)
plt.plot(radiuslist,'.')
plt.title('orbital separation iterative solution')
plt.ylabel('r')
plt.xlabel('Iteration')
plt.show()

#r the vector and its time derivative

asteroid = rholist[len(rholist)-1]*s2
R1 = earth[0]
R2 = earth[1]
R3 = earth[2]

r = R2 + asteroid

rmag = radiuslist[len(radiuslist)-1]

R2dot = ((tau3/(tau1*(tau1+tau3)))*(R2-R1)) + ((tau1/(tau3*(tau1+tau3)))*(R3-R2))

rdot = R2dot + rholist[len(rholist)-1]*s2dot + rhovel*s2

V = np.linalg.norm(rdot)

a = (rmag*k**2)/((2*k**2)-(rmag*V**2))

h = np.cross(r,rdot)

hx,hy,hz = h

hmag = np.linalg.norm(h)

Omega = np.arctan(-hx/hy)
i = np.arccos(hz/hmag)
e = np.sqrt(1-((hmag**2)/(a*k**2)))
E = np.arccos((a-rmag)/(a*e))
E = -E
nu = trueanomaly(e,E)
M = meananomalyE(E,e)
sumnu = np.arccos((r[0]*np.cos(Omega)+r[1]*np.sin(Omega))/rmag)
omega = sumnu-nu
n = np.sqrt((k**2)/(a**3))
tau = Julian[1]-(M/n)