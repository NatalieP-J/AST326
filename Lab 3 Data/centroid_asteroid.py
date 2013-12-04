import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

asteroid_x = [1089.5,1090.5,1065,1088.5,1071.5]
asteroid_y = [977.5,1039.5,1008,1039,1025]

centroids = True

skies = []
darks = []
flats = []
corrected = []

days = [19,20,22,23,28]

thresholds = []
backgrounds = []

actual_ax = []
actual_ay = []


for urania in range(4):
	print 'Working on {0}'.format(urania)
	skyname = '30Urania-S001-R001-C001-r-{0}.fts'.format(urania)
	darkname = 'Dark-S001-R001-C001-B2-{0}.fts'.format(urania)
	flatname = 'AutoFlat-Dusk-r-Bin2-001-{0}.fts'.format(urania)
	s = pf.open(skyname)
	sky_time = s[0].header['exptime']
	sky_date = s[0].header['date']
	print 'Raw exposure {0}s'.format(sky_time)
	print 'Raw date {0}'.format(sky_date)
	sky = pf.getdata(skyname)
	sky /= sky_time
	s = pf.open(darkname)
	dark_time = s[0].header['exptime']
	print 'Dark exposure {0}s'.format(dark_time)
	dark = pf.getdata(darkname)
	dark /= dark_time
	s = pf.open(flatname)
	flat_time = s[0].header['exptime']
	print 'Flat exposure {0}s'.format(flat_time)
	flat = pf.getdata(flatname)
	flat /= flat_time
	skies.append(sky)
	darks.append(dark)
	flats.append(flat)
	f = flat-dark
	f = f/np.median(f)
	correct = ((sky-dark)/f)[::-1]
	corrected.append(correct)

	if centroids == True:

		print 'Finding local maxima'
		
		max_x = []
		max_y = []
		intensities = np.zeros((2048,2048))
		threshold = 100
		for i in range(1,len(correct)-1):
			row = correct[i]
			for j in range(1,len(row)-1):
				if row[j] >= threshold:
					if row[j] >= row[j-1] and row[j] >= row[j+1] and row[j] >= correct[i-1][j] and row[j] >= correct[i+1][j] :
						max_x.append(j)
						max_y.append(i)
						intensities[i][j] = row[j]


		#FIND THE CENTROIDS
		background = 50
		radius = 25

		print 'Finding centroids'

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

		np.savetxt('Jan{0}_Xcentroid.txt'.format(days[urania]),centroid_x)
		np.savetxt('Jan{0}_Ycentroid.txt'.format(days[urania]),centroid_y)

		a_x = asteroid_x[urania]
		a_y = asteroid_y[urania]

		separation = []
		for i in range(len(centroid_x)):
			xdiff = centroid_x[i]-a_x
			ydiff = centroid_y[i]-a_y
			R = (xdiff**2)+(ydiff**2)
			separation.append(R)

		for j in range(len(separation)):
			if separation[j] == min(separation):
				asteroid = j
				actual_ax.append(centroid_x[j])
				actual_ay.append(centroid_y[j])

		plt.subplot(2,3,urania+1)
		plt.imshow(correct,cmap='gray_r',origin = 'lower',norm = LogNorm(vmin = 10,vmax = 1000))
		plt.colorbar()
		plt.plot(a_x,a_y,'g.')
		plt.plot(centroid_x[asteroid],centroid_y[asteroid],'hb')
		plt.xlim(0,2048)
		plt.ylim(0,2048)
		for i in range(len(centroid_x)):
			if centroid_x[i] != centroid_x[asteroid]:
			    circ  = plt.Circle((centroid_x[i],centroid_y[i]),radius = 18, color='r',fill=False)
			    plt.gcf()
			    plt.gca().add_artist(circ)
		plt.xlabel('x [pixel]')
		plt.ylabel('y [pixel]')
		plt.title('30 Urania January {0}'.format(days[urania]))

actual_ay.pop(2)
actual_ax.pop(2)

centroids4 = True

urania = 4
print 'Working on {0}'.format(urania)
skyname = '30Urania-S001-R001-C001-r-{0}.fts'.format(urania)
darkname = 'Dark-S001-R001-C001-B2-{0}.fts'.format(urania)
flatname = 'AutoFlat-Dusk-r-Bin2-001-{0}.fts'.format(urania)
s = pf.open(skyname)
sky_time = s[0].header['exptime']
sky_date = s[0].header['date']
print 'Raw exposure {0}s'.format(sky_time)
print 'Raw date {0}'.format(sky_date)
sky = pf.getdata(skyname)
sky /= sky_time
s = pf.open(darkname)
dark_time = s[0].header['exptime']
print 'Dark exposure {0}s'.format(dark_time)
dark = pf.getdata(darkname)
dark /= dark_time
s = pf.open(flatname)
flat_time = s[0].header['exptime']
print 'Flat exposure {0}s'.format(flat_time)
flat = pf.getdata(flatname)
flat /= flat_time
skies.append(sky)
darks.append(dark)
flats.append(flat)
f = flat-dark
f = f/np.median(f)
correct = ((sky-dark)/f)[::-1]
corrected.append(correct)

if centroids4 == True:

	print 'Finding local maxima'
	
	max_x = []
	max_y = []
	intensities = np.zeros((2048,2048))
	threshold = 600
	for i in range(1,len(correct)-1):
		row = correct[i]
		for j in range(1,len(row)-1):
			if row[j] >= threshold and j < 2000:
				if row[j] >= row[j-1] and row[j] >= row[j+1] and row[j] >= correct[i-1][j] and row[j] >= correct[i+1][j] :
					max_x.append(j)
					max_y.append(i)
					intensities[i][j] = row[j]


	#FIND THE CENTROIDS
	background = 300
	radius = 25

	print 'Finding centroids'

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

	np.savetxt('Jan{0}_Xcentroid.txt'.format(days[urania]),centroid_x)
	np.savetxt('Jan{0}_Ycentroid.txt'.format(days[urania]),centroid_y)

	a_x = asteroid_x[urania]
	a_y = asteroid_y[urania]

	separation = []
	for i in range(len(centroid_x)):
		xdiff = centroid_x[i]-a_x
		ydiff = centroid_y[i]-a_y
		R = (xdiff**2)+(ydiff**2)
		separation.append(R)

	for j in range(len(separation)):
		if separation[j] == min(separation):
			asteroid = j
			actual_ax.append(centroid_x[j])
			actual_ay.append(centroid_y[j])

	plt.subplot(2,3,urania+1)
	plt.imshow(correct,cmap='gray_r',origin = 'lower',norm = LogNorm(vmin = 300,vmax = 1000))
	plt.colorbar()
	plt.plot(a_x,a_y,'g.')
	plt.plot(centroid_x[asteroid],centroid_y[asteroid],'hb')	
	plt.xlim(0,2048)
	plt.ylim(0,2048)
	for i in range(len(centroid_x)):
		if centroid_x[i] != centroid_x[asteroid]:
		    circ  = plt.Circle((centroid_x[i],centroid_y[i]),radius = 18, color='r',fill=False)
		    plt.gcf()
		    plt.gca().add_artist(circ)
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	plt.title('30 Urania January {0}'.format(days[urania]))
	plt.show()

