import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

centroids = False


#Import NGC7331 Data

galaxy = 'NGC7331-S001-R001-C001-r.fts' #exposure time 240s
dark = 'Dark-S001-R003-C003-B2.fts'
flat = 'combflatr.fits'
d = pf.getdata(dark)
x = pf.getdata(galaxy)
f = pf.getdata(flat)
fnew = f/np.median(f) #get appropriate variation in flat field
correct = ((x-d)/fnew)[::-1] #reduce the raw data

plt.figure()
plt.imshow(fnew,origin='lower',cmap = 'winter',vmin = 0.7,vmax=1.2)
plt.title('Normalized Flat Spectrum for NGC7331')
plt.xlabel('x [pixel]')
plt.ylabel('y [pixel]')
plt.colorbar()
plt.show()

#zoom in on a relatively uniform section of the flat and plot
#pixel to pixel variations
fycrop = fnew[1350:1360]
fzoom = np.zeros([10,10])
for i in range(len(fycrop)):
	for j in range(630,640):
		fzoom[i][j-630] = fycrop[i][j]

plt.figure()
plt.imshow(fzoom,origin='lower',cmap = 'winter',vmin = 1.07,vmax=1.12)
plt.title('Normalized Flat Spectrum for NGC7331')
plt.xlabel('x [pixel-630]')
plt.ylabel('y [pixel-1350]')
plt.colorbar()
plt.show()

mean = np.mean(correct)
sigma = np.std(correct)

if centroids == True:

	#find local maxima
	#leave out the edge values, as we will not be able to get good centroids for them
	#thus they will not be useful

	max_x = [] #array to hold xpixel locations
	max_y = [] #array to hold ypixel locations

	intensities = np.zeros((2048,2048)) #array to hold intensities

	threshold = 35000 #lower limit on what is considered a star

	for i in range(1,len(correct)-1): #ignore the top and bottom rows
		#extract a row of the field
		row = correct[i] 
		
		for j in range(1,len(row)-1): #check through all values in the row except the edges
			if row[j] >= threshold: #if the intensity is greater than our threshold
				#and is greater than or equal to the intesity of the four surrounding pixels
				if row[j] >= row[j-1] and row[j] >= row[j+1] and row[j] >= correct[i-1][j] and row[j] >= correct[i+1][j] :
					max_x.append(j) #add its x location
					max_y.append(i) #add its y location
					intensities[i][j] = row[j]


	#FIND THE CENTROIDS
	background = mean+2*sigma #background limit
	radius = 25 #pixel radius around each star, in which no other stars seems to intrude

	centroid_x = [] #array to hold x centroids
	centroid_y = [] #array to hold y centroids

	#begin checking through each of the located loca maxima
	for k in range(len(max_x)):
		peak = correct[max_y[k]][max_x[k]] #intensity location at local max
		x_sum = 0 #x centroid initial value
		y_sum = 0 #y centroid initial value
		Ix_sum = 1 #intensity intial value - set to 1 because I was having divide by zero problems
		Iy_sum = 1 #since intensities are of order at least 10000, 1 doesn't make much difference
		j = -25 
		#check a square with radius 25 around our local mix
		while -25<=j<=25:
			n=-25
			while -25<=n<=25:
				try:
				 	xlim = correct[max_y[k]+n][max_x[k]+j]
				 	ylim = correct[max_y[k]+j][max_x[k]+n]

				 	#if the pixel exceeds background value, add it to the sum
					if correct[max_y[k]+n][max_x[k]+j] > background:
						x_sum += correct[max_y[k]+n][max_x[k]+j]*(max_x[k]+j)
						Ix_sum += correct[max_y[k]+n][max_x[k]+j]
					if correct[max_y[k]+j][max_x[k]+n] > background:
						y_sum += correct[max_y[k]+j][max_x[k]+n]*(max_y[k]+j)
						Iy_sum += correct[max_y[k]+j][max_x[k]+n]
					n+=1

				#if a pixel is not available with in that radius, move on to the next one
				except IndexError:
					j=26
					n=26
			j+=1

		#construct the weighted averages
		centroid_x.append(x_sum/Ix_sum)
		centroid_y.append(y_sum/Iy_sum)

	#save the centroid values
	np.savetxt('Xcentroid.txt',centroid_x)
	np.savetxt('Ycentroid.txt',centroid_y)

if centroids == False:

	centroid_x = np.loadtxt('Xcentroid.txt')
	centroid_y = np.loadtxt('Ycentroid.txt')

#use estimated location of galaxy to find it in the centroid list
for i in range(len(centroid_x)):
	if abs(1178.32-centroid_x[i])<5:
		index = i

#manually remove galaxy
centroid_x = np.delete(centroid_x,index)
centroid_y = np.delete(centroid_y,index)

#plot the galaxy image
plt.figure()
plt.imshow(correct,cmap='gray_r',origin = 'lower',norm = LogNorm(vmin = 3000,vmax = 100000))
plt.colorbar()
for i in range(len(centroid_x)):
    circ  = plt.Circle((centroid_x[i],centroid_y[i]),radius = 18, color='r',fill=False)
    plt.gcf()
    plt.gca().add_artist(circ)
plt.xlabel('x [pixel]')
plt.ylabel('y [pixel]')
plt.title('Corrected CCD Image of Field Surrounding NGC7331')
plt.xlim(0,2048)
plt.ylim(0,2048)
plt.show()

