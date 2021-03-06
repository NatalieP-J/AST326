import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

array = np.array

rad = np.pi/180.

k = 0.01720209895

expected_r_Jan23=[5.557860055930863E-01,2.046029478081670E+00,6.190352100780391E-02]

#FUNCTIONS 

def rho(k,R,s,sdot,sdouble,r):
	numerator = np.dot(sdot,np.cross(R,s))
	denominator = np.dot(sdot,np.cross(sdouble,s))
	Rmag = np.linalg.norm(R)
	constants = (k**2)*((1./Rmag**3)-(1./r**3))
	return constants*(numerator/denominator)

def rhoerr(k,R,s,sdot,sdouble,serr,sdoterr,sdoubleerr,r,rerr):
	Rmag = np.linalg.norm(R)
	C = ((k**2)/2)*((1./Rmag**3)-(1./r**3))
	deriv = np.array([1,1,1])
	delrhos = C*((np.dot(sdot,np.cross(R,deriv))*np.dot(sdot,np.cross(sdouble,s))-np.dot(sdot,np.cross(sdouble,deriv))*np.dot(sdot,np.cross(R,s)))/(np.dot(sdot,np.cross(sdouble,s)))**2)
	delrhosdot = C*((np.dot(deriv,np.cross(R,s)*np.dot(sdot,np.cross(sdouble,s))-np.dot(deriv,np.cross(sdouble,s))*np.dot(sdot,np.cross(R,s))))/(np.dot(sdot,np.cross(sdouble,s)))**2)
	delrhosdou = C*((-np.dot(sdot,np.cross(R,s))*np.dot(sdot,np.cross(deriv,s)))/(np.dot(sdot,np.cross(sdouble,s)))**2)
	delrhor = (k**2)*(3/r**4)*((np.dot(sdot,np.cross(R,s)))/np.dot(sdot,np.cross(sdouble,s)))
	return np.sqrt((delrhos*serr)**2+(delrhosdot*sdoterr)**2+(delrhosdou*sdoubleerr)**2+(delrhor*rerr))

def rhodot(k,R,s,sdot,sdouble,r):
	numerator = np.dot(sdouble,np.cross(R,s))
	denominator = np.dot(sdouble,np.cross(sdot,s))
	Rmag = np.linalg.norm(R)
	constants = ((k**2)/2)*((1./Rmag**3)-(1./r**3))
	return constants*(numerator/denominator)

def rhodoterr(k,R,s,sdot,sdouble,serr,sdoterr,sdoubleerr,r,rerr):
	Rmag = np.linalg.norm(R)
	C = ((k**2)/2)*((1./Rmag**3)-(1./r**3))
	deriv = np.array([1,1,1])
	delrhodotr = ((3*k**2)/(2*r**4))*(np.dot(sdouble,np.cross(R,s))/np.dot(sdouble,np.cross(sdot,s)))
	delrhodots = C*((np.dot(sdouble,np.cross(R,deriv))*np.dot(sdouble,np.cross(sdot,s))-np.dot(sdouble,np.cross(R,s))*np.dot(sdouble,np.cross(sdot,deriv)))/((np.dot(sdouble,np.cross(sdot,s)))**2))
	delrhodotsdot = C*((-np.dot(sdouble,np.cross(R,s))*np.dot(sdouble,np.cross(deriv,s)))/((np.dot(sdouble,np.cross(sdot,s)))**2))
	delrhodotsdou = C*((np.dot(deriv,np.cross(R,s))*np.dot(sdouble,np.cross(sdot,s))-np.dot(deriv,np.cross(sdot,s))*np.dot(sdouble,np.cross(R,s)))/((np.dot(sdouble,np.cross(sdot,s)))**2))
	return np.sqrt((delrhodotr*rerr)**2+(delrhodots*serr)**2+(delrhodotsdot*sdoterr)**2+(delrhodotsdou*sdoubleerr)**2)

def radius(rhoval,R,s):
	Rmag = np.linalg.norm(R)
	return np.sqrt((rhoval**2)+(Rmag**2)+2*rhoval*np.dot(R,s))

def radiuserr(rhoval,rhovalerr,R,s,serr):
	Rmag = np.linalg.norm(R)
	r = np.sqrt((rhoval**2)+(Rmag**2)+2*rhoval*np.dot(R,s))
	delrrho = (1./r)*(rhoval+np.dot(R,s))
	delrs = (1./r)*(rhoval*np.dot(R,s))
	return np.sqrt((delrrho*rhovalerr)**2+(delrs*serr)**2)

def trueanomaly(eccentricity,angles):
	return 2*np.arctan(np.sqrt((1+eccentricity)/(1-eccentricity))*np.tan(angles/2))

def trueanomaly1(eccentricity,angles):
	return 2*np.arctan(np.sqrt(-(1+eccentricity)/(1-eccentricity))*np.tan(angles/2))

def meananomalyE(E,e):
	return E-e*np.sin(E)

def ellipseradius(angles,semimajor,eccentricity):
	return semimajor*(1.0-eccentricity*np.cos(angles))

def meananomalytime(n,tau,t):
	return n*(t-tau)

def circular(angles,radius):
	x = radius*np.cos(angles)
	y = radius*np.sin(angles)
	return x,y

def rvecerr(R,rhoval,s,rhoerr,serr):
	return np.sqrt((s*rhoerr)**2+(rhoval*serr)**2)

def rdotvecerr(R,rhoval,s,rhovel,sdot,rhoerr,serr,rhovelerr,sdoterr):
	return np.sqrt((sdot*rhoerr)**2+(rhoval*sdoterr)**2+(rhovel*serr)**2+(s*rhovelerr)**2)

#FIXED PLATE CONSTANTS 4

alpha = array([ 44.45521841,  44.68467768,  45.17098658,  45.40514157,  46.7449028 ])
delta = array([ 19.21813581,  19.27914174,  19.36052613,  19.40083071,  19.62835052])

alpha *= rad
delta *= rad

#ERRORS IN MEASUREMENTS

alphaerr = array([ 0.00577571,  0.0058012 ,  0.00582736,  0.00571606,  0.00588477])
deltaerr = array([ 0.0059024 ,  0.0055641 ,  0.00571975,  0.00558783,  0.00578666])

xeq = np.cos(alpha)*np.cos(delta)
yeq = np.sin(alpha)*np.cos(delta)
zeq = np.sin(delta)

def coordshifterr(alpha,delta,alphaerr,deltaerr):
	xeq = np.cos(alpha)*np.cos(delta)
	yeq = np.sin(alpha)*np.cos(delta)
	zeq = np.sin(delta)
	xerr = np.sqrt((-yeq*alphaerr)**2+(-np.sin(delta)*np.cos(alpha)*deltaerr)**2)
	yerr = np.sqrt((xeq*alphaerr)**2+(-np.sin(alpha)*np.sin(delta)*deltaerr)**2)
	zerr = np.sqrt((np.cos(delta)*deltaerr)**2)
	return xerr,yerr,zerr

xerr,yerr,zerr = coordshifterr(alpha,delta,alphaerr,deltaerr)

epsilon = 23.43929111*rad

T = [[1,0,0],[0,np.cos(epsilon),np.sin(epsilon)],[0,-np.sin(epsilon),np.cos(epsilon)]]
T = np.matrix(T)

s_err = []
s = []

for i in range(len(alpha)):
	r = [xeq[i],yeq[i],zeq[i]]
	r_err = [xerr[i],yerr[i],zerr[i]]
	r = np.matrix(r)
	r_err = np.matrix(r_err)
	r = r.T
	r_err = r_err.T
	ecliptic = T*r
	ecliptic_err = T*r_err
	ecliptic = np.array(ecliptic.T)
	ecliptic_err = np.array(ecliptic_err.T)
	s.append([ecliptic[0][0],ecliptic[0][1],ecliptic[0][2]])
	s_err.append([ecliptic_err[0][0],ecliptic_err[0][1],ecliptic_err[0][2]])

s = np.array(s)
s_err = np.array(s_err)

earth = np.array([	[-4.852325900266916e-01,  8.563819442309909e-01, -2.676783649804444e-05],
					[-5.005639235484095e-01,  8.476786115155864e-01, -2.723171022212407e-05],
					[-5.311509164664779e-01,  8.292136310803589e-01, -2.794675828254335e-05],
					[-5.450861676152123e-01,  8.202934515047903e-01, -2.812109168370201e-05],
					[-6.143806933050102e-01,  7.707734163615326e-01, -2.700019419225799e-05]])


actuals = np.array([[1.073459020670066E+00,  1.175275572464343E+00,  6.266997776117721E-02],
					[1.077052532037185E+00,  1.188320067084311E+00,  6.242712185034433E-02],
					[1.083793604385678E+00,  1.215389298506962E+00,  6.192856377902563E-02],
					[1.086660494660806E+00,  1.228204989345288E+00,  6.169477872353334E-02],
					[1.098706174782393E+00,  1.296878129932550E+00,  6.046177329415765E-02]])

snew = []
rhovals = []

for i in range(len(actuals)):
	snew.append(actuals[i]/np.linalg.norm(actuals[i]))
	rhovals.append(np.linalg.norm(actuals[i]))

times = [	'2012-01-20 04:28:30',
			'2012-01-21 04:40:27',
			'2012-01-23 05:43:40',
			'2012-01-24 04:26:48',
			'2012-01-29 01:27:18']
t = Time(times,format='iso',scale = 'utc')

Julian = t.jd

threesome = [1,2,3,4]

s1 = np.array(s[threesome[0]])
s2 = np.array(s[threesome[1]])
s3 = np.array(s[threesome[2]])
s1err = np.array(s_err[threesome[0]])
s2err = np.array(s_err[threesome[1]])
s3err = np.array(s_err[threesome[2]])

tau1 = Julian[threesome[1]]-Julian[threesome[0]]
tau3 = Julian[threesome[2]]-Julian[threesome[1]]

s2dot = ((tau3/(tau1*(tau1+tau3)))*(s2-s1)) + ((tau1/(tau3*(tau1+tau3)))*(s3-s2))
s2doterr = ((tau3/(tau1*(tau1+tau3)))*(s2err-s1err)) + ((tau1/(tau3*(tau1+tau3)))*(s3err-s2err))
s2double = ((2/(tau3*(tau1+tau3)))*(s3-s2)) - ((2/(tau1*(tau1+tau3)))*(s2-s1))
s2doubleerr = ((2/(tau3*(tau1+tau3)))*(s3err-s2err)) - ((2/(tau1*(tau1+tau3)))*(s2err-s1err))

#s2doterr = np.sqrt(((tau3/(tau1*(tau1+tau3)))*s1err)**2+((tau1/(tau3*(tau1+tau3)))*s3err)**2+(((tau3/(tau1*(tau1+tau3)))-(tau1/(tau3*(tau1+tau3))))*s2err)**2)
#s2doubleerr = np.sqrt(((2/(tau1*(tau1+tau3)))*s1err)**2+((2/(tau3*(tau1+tau3)))*s3err)**2+(((-2/(tau3*(tau1+tau3)))+(-2/(tau1*(tau1+tau3))))*s2err)**2)

mags2err = np.linalg.norm(s2err)
mags2doterr = np.linalg.norm(s2doterr)
mags2doubleerr = np.linalg.norm(s2doubleerr)

r0 = 2.0 #AU
radiuslist = []
radiuslisterr = []
rholist = []
rholisterr = []
radiuslist.append(r0)
radiuslisterr.append(0)

R2 = earth[threesome[1]]

for i in range(100):
	p = rho(k,R2,s2,s2dot,s2double,radiuslist[i])
	rholist.append(p)
	perr = rhoerr(k,R2,s2,s2dot,s2double,mags2err,mags2doterr,mags2doubleerr,radiuslist[i],radiuslisterr[i])
	rholisterr.append(perr)
	rval = radius(rholist[i],R2,s2)
	radiuslist.append(rval)
	rvalerr = radiuserr(rholist[i],rholisterr[i],R2,s2,mags2err)
	radiuslisterr.append(rvalerr)

rhovel = rhodot(k,R2,s2,s2dot,s2double,radiuslist[len(rholist)-1])
rhovelerr = rhodoterr(k,R2,s2,s2dot,s2double,mags2err,mags2doterr,mags2doubleerr,radiuslist[len(rholist)-1],radiuslisterr[len(rholist)-1])

plt.figure()
plt.subplot(211)
plt.plot(rholist,'.')
plt.title(r'$\rho$ iterative solution')
plt.ylabel(r'$\rho$')
plt.xlabel('Iteration')
plt.subplot(212)
plt.plot(radiuslist,'.')
plt.title('radius iterative solution')
plt.ylabel('r')
plt.xlabel('Iteration')
#plt.show()

asteroid = rholist[len(rholist)-1]*s2
R1 = earth[threesome[0]]
R2 = earth[threesome[1]]
R3 = earth[threesome[2]]

deriv = [1,1,1]

r = R2 + asteroid
rerr = rvecerr(R2,rholist[len(rholist)-1],s2,rholisterr[len(rholist)-1],s2err)

rmag = radiuslist[len(radiuslist)-1]

R2dot = ((tau3/(tau1*(tau1+tau3)))*(R2-R1)) + ((tau1/(tau3*(tau1+tau3)))*(R3-R2))

rdot = R2dot + rholist[len(rholist)-1]*s2dot + rhovel*s2
rdoterr = rdotvecerr(R2,rholist[len(rholist)-1],s2,rhovel,s2dot,rholisterr[len(rholist)-1],s2err,rhovelerr,s2doterr)

V = np.linalg.norm(rdot)
Verr = np.linalg.norm(rdoterr)

a = (rmag*k**2)/((2*k**2)-(rmag*V**2))
aerr = np.sqrt((((k**2)*(rmag**2)*V)/(2*k**2-rmag*V**2))**Verr)

h = np.cross(r,rdot)
herr = np.sqrt((np.cross(deriv,rdot)*rerr)**2+(np.cross(r,deriv)*rdoterr)**2)

hx,hy,hz = h
hxerr,hyerr,hzerr = herr

hmag = np.linalg.norm(h)
hmagerr = np.linalg.norm(herr)

Omega = np.arctan2(-hx,hy)
Omegaerr = np.arctan2(-hxerr,hyerr)
Omega += np.pi
inc = np.arccos(hz/hmag)
incerr = np.arccos(hzerr/hmagerr)
e = np.sqrt(1-((hmag**2)/(a*k**2)))
e_err = np.sqrt(-(1-((hmagerr**2)/(aerr*k**2))))
E = np.arccos((a-rmag)/(a*e))
Eerr = np.arccos((aerr-radiuslisterr[9])/(aerr*e_err))
Eval = E
nu = trueanomaly(e,E)
nuerr = trueanomaly1(e_err,Eerr)
M = meananomalyE(E,e)
Merr = meananomalyE(Eerr,e_err)
sumnu = np.arccos((r[0]*np.cos(Omega)+r[1]*np.sin(Omega))/rmag)
sumnuerr = np.arccos((rerr[0]*np.cos(Omegaerr)+rerr[1]*np.sin(Omegaerr))/radiuslisterr[9])
omega = sumnu-nu
n = np.sqrt((k**2)/(a**3))
tau = Julian[threesome[1]]-(M/n)
tauerr = (Merr/n)

period = (2*np.pi)/n

#Generate time intervals
start = tau
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

theta = v+omega

#FIGURE 3

plt.figure()
plt.plot(J,E,'k',label='Eccentric Anomaly')
plt.plot(J,M,'k:',label='Mean Anomaly')
plt.xlabel('Julian Day - 2450000')
plt.ylabel('Anomaly [radians]')
plt.title('Anomaly Over Time')
plt.legend(loc='best')
#plt.show()

#FIGURE 4

plt.figure()
plt.subplot(211)
plt.plot(J,r,'k')
plt.xlabel('Julian Day - 2450000')
plt.ylabel('r [AU]')
plt.title('Orbital Separation Over Time')
plt.subplot(212)
plt.plot(J,v,'k')
plt.xlabel('Julian Day - 2450000')
plt.ylabel('v [radians]')
plt.title('True Anomaly Over Time')
plt.tight_layout()
#plt.show()

#FIGURE 5

plt.figure()
plt.plot(circular(theta,a)[0],circular(theta,a)[1],':',label = 'Circular orbit radius = {0}'.format(np.round(a,2)))
plt.plot(r*np.cos(theta),r*np.sin(theta),'k',label = 'Calculated orbit')
plt.plot(0,0,'k+')
plt.plot(r[0]*np.cos(theta[0]),r[0]*np.sin(theta[0]),'ko')
plt.axis('equal')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.legend()
plt.title('Position of Urania in its Orbital Plane')
#plt.show()

time = Julian[threesome[3]]

Rearth = earth[threesome[3]]

Rearth = [-6.172959176680726E-01,  7.674840479224522E-01, -2.337958345652992E-05]

for i in range(len(JulianDay)):
	if np.round(JulianDay[i])==np.round(time):
		index = i

Mnew = M[index]

Enew = E[index]

v = trueanomaly(e,Enew)

theta = v+omega
thetaerr = np.sqrt(nuerr**2 + sumnuerr**2)

rmag = a*(1-e*np.cos(Enew))
rmagerr = aerr*(1-e_err*np.cos(Eerr))


rvec = [rmag*np.cos(theta),rmag*np.sin(theta),0]
rvecerr = [rmagerr*np.cos(thetaerr),rmagerr*np.sin(thetaerr),0]

TzTx = [[np.cos(Omega),-np.sin(Omega)*np.cos(inc),np.sin(Omega)*np.sin(inc)],
		[np.sin(Omega),np.cos(Omega)*np.cos(inc),-np.cos(Omega)*np.sin(inc)],
		[0,np.sin(inc),np.cos(inc)]]

TzTxerr = [[np.cos(Omegaerr),-np.sin(Omegaerr)*np.cos(incerr),np.sin(Omegaerr)*np.sin(incerr)],
		[np.sin(Omegaerr),np.cos(Omegaerr)*np.cos(incerr),-np.cos(Omegaerr)*np.sin(incerr)],
		[0,np.sin(incerr),np.cos(incerr)]]

TzTx = np.matrix(TzTx)
TzTxerr = np.matrix(TzTxerr)

rvec = np.matrix(rvec)
rvecerr = np.matrix(rvecerr)

r_ecliptic = TzTx*rvec.T
r_ecliptic_err = TzTxerr*rvecerr.T

rhos = np.array(r_ecliptic.T) - Rearth
rhoserr = np.array(r_ecliptic_err.T)

rhos = T.T*np.matrix(rhos).T

rhos = np.array(rhos.T)

rhoserr = np.array((T.T*np.matrix(rhoserr).T).T)

rho = np.linalg.norm(rhos)

s = rhos/rho

x,y,z = rhos[0]
xerr,yerr,zerr = rhoserr[0]

alpha = np.arctan2(y,x)
alphaerr = np.arctan2(yerr,xerr)

delta = np.arcsin(z/rho)
deltaerr = np.arcsin(zerr/np.linalg.norm(rhoserr))
