import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np

#350936

galaxy = 'NGC7331-S001-R001-C001-r.fts' #exposure time 240s
dark = 'Dark-S001-R003-C003-B2.fts'
flat = 'combflatr.fits'
d = pf.getdata(dark)
dhdr = pf.getheader(dark)
x = pf.getdata(galaxy)
f = pf.getdata(flat)
fnew = f/np.median(f)
correct = ((x-d)/fnew)[::-1]
hdr = pf.getheader(galaxy)
hdulist = pf.open(galaxy)

max_x = []
max_y = []
intensities = np.zeros((2048,2048))
threshold = 35000
for i in range(1,len(correct)-1):
	row = correct[i]
	for j in range(1,len(row)-1):
		if row[j] >= threshold:
			if row[j] >= row[j-1] and row[j] >= row[j+1] and row[j] >= correct[i-1][j] and row[j] >= correct[i+1][j] :
				max_x.append(j)
				max_y.append(i)
				intensities[i][j] = row[j]


#FIND THE CENTROIDS
background = 10000
radius = 25

centroid_x = []
centroid_y = []

for k in range(len(max_x)):
	print 'Working on peak {0} of {1}'.format(k,len(max_x))
	peak = correct[max_y[k]][max_x[k]]
	x_sum = 0
	y_sum = 0
	Ix_sum = 1
	Iy_sum = 1
	j = -25
	while -25<=j<=25:
		n=-25
		while -25<=n<=25:
			try:
			 	xlim = correct[max_y[k]+n][max_x[k]+j]
			 	ylim = correct[max_y[k]+j][max_x[k]+n]
				if correct[max_y[k]+n][max_x[k]+j] > background:
					x_sum += correct[max_y[k]+n][max_x[k]+j]*(max_x[k]+j)
					Ix_sum += correct[max_y[k]+n][max_x[k]+j]
				if correct[max_y[k]+j][max_x[k]+n] > background:
					y_sum += correct[max_y[k]+j][max_x[k]+n]*(max_y[k]+j)
					Iy_sum += correct[max_y[k]+j][max_x[k]+n]
				n+=1
			except IndexError:
				j=26
				n=26
		j+=1
	centroid_x.append(x_sum/Ix_sum)
	centroid_y.append(y_sum/Iy_sum)

np.savetxt('Xcentroid.txt',centroid_x)
np.savetxt('Ycentroid.txt',centroid_y)

plt.figure()
plt.imshow(correct,cmap='gray_r',origin = 'lower',norm = LogNorm(vmin = 10,vmax = 100000))
plt.colorbar()
plt.plot(centroid_x,centroid_y,'.')
plt.xlim(0,2048)
plt.ylim(0,2048)
#plt.scatter(centroid_x,centroid_y,s=50,facecolors='none',edgecolors='r')
plt.show()
