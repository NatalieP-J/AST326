import matplotlib.pyplot as plt
import numpy as np
from astropy.stats.funcs import sigma_clip

files = '/h/ungrad1/stostad/lab2/neon/neon000'

def sigmas(files, nofiles = 48, plot = False):
    '''Returns a mean spectrum in the form of a 2d array - one dimension
    is the wavelengths, the other is the flux files.

    The input is the start of the filenames; i.e. if you have 1040 data
    files named "mortenisthebest00001.txt", "mortenisthebest00002.txt" etc.,
    the input is files="mortenisthebest0" (the file format is assumed to
    be .txt) and nofiles = 1040.
    '''
    all_data = np.zeros((nofiles, 2048, 2))
    # First loop over the files to get data from each and put the data
    # in all_data:
    for indx in range(nofiles):
        file1 = files + str(indx).zfill(len(str(nofiles))) + '.txt'
        data1 = np.genfromtxt(file1, skip_header = 17, skip_footer = 1)
        all_data[indx] = data1
    datamean = np.zeros((2,2048))
    # We want to have a mean array datamean that gives us the x-values
    # and mean flux values for each wavelength in the spectrum.

    # The x-values are assumed to be the same for every spectrum.
    datamean[0] = all_data[0,:,0]

    # Do a sigma clip to take away outliers
    for i in range(len(all_data[0,:,0])):
        bool_clip = sigma_clip(all_data[:,i,1], sig = 1.5)[1]
        datamean[1,i] = np.mean(all_data[:,i,1][bool_clip])
    # datamean[1] = np.mean(all_data[:,:,1], axis = 0)
    if plot == True:
        plt.clf()
        plt.plot(datamean[0], datamean[1])
        plt.show()
        plt.semilogy()
        plt.axis([530,680,15,60000])
        plt.xlabel('Wavelength (nm)', fontsize = 20)
        plt.ylabel('Flux (counts)', fontsize = 20)
        plt.title('The raw mean of the spectra', fontsize = 25)
    return datamean

def maxim(file1, threshold = 800):
    '''Finds the maxima of a spectrum (which is given in the format of
    the data returned from sigmas()).
    '''
    if type(file1) == str:
        data = np.genfromtxt(file1, skip_header = 17, skip_footer = 1)
    else:
        data = file1
    maxima = None
    for i in range(len(data[0,1:-1])):
        if data[1,i] > data[1,i-1] and data[1,i] > data[1,i+1] and data[1,i] > threshold:
            if maxima == None:
                maxima = np.array(data[0,i])
            else:
                maxima = np.vstack((maxima, np.array(data[0,i])))
            print i, data[:,i]
    np.savetxt('maxima.dat', maxima)

def maxima(files):
    '''Uses both functions to find the sigma clipped mean and then find
    the maxima of that mean spectrum - all from just a list of files.
    '''
    mean = sigmas(files)
    maxim(mean)
