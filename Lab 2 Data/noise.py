import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq

files = 'October25/lamp_Oct25_3_0'

dark_files = 'October25/dark_Oct25_0'

def guess(params,independent):
    return params[0]*independent+params[1]

def residuals(params,independent,theory):
    return abs(theory - guess(params,independent))

def polyguess(params,independent):
	return params[0]*independent+params[1]*independent**2+params[2]*independent**3+params[3]

def polyresiduals(params,independent,theory):
	return abs(theory - polyguess(params,independent))

def import_files(files, nofiles = 1024):
    all_data = np.zeros((nofiles, 2048, 2))
    for indx in range(nofiles):
        file1 = files + str(indx).zfill(len(str(nofiles))) + '.txt'
        data1 = np.genfromtxt(file1, skip_header = 17, skip_footer = 1)
        all_data[indx] = data1
    return all_data

def plot_pixel(all_data,ra):
    for j in range(len(ra)):
        plt.subplot(len(ra)/2,len(ra)/2,j+1)
        pixel = []
        for i in range(len(all_data)):
            pixel.append(all_data[i][:,1][ra[j]])
        x = range(0,len(pixel))
        plt.xlabel('Sample Number')
        plt.ylabel('Intensity (ADU)')
        plt.title('Pixel Number {0}'.format(ra[j]))
        plt.plot(x,pixel)
        plt.xlim(0,1024)
    plt.tight_layout()
    plt.show()

def noise(all_data):
    all_data = all_data[800:1000]

    means = np.mean(all_data[:,:,1], axis = 0)
    stds = np.std(all_data[:,:,1], axis = 0)

    means_lsq = means[np.where(means < 40000)]
    stds_lsq = stds[np.where(means < 40000)]

    means_fit = means_lsq[np.where(means_lsq > 100)]
    stds_fit = stds_lsq[np.where(means_lsq > 100)]

    params0 = [0.2,1,0.1,40]

    p, pcov_x, info, mesg, success = leastsq(polyresiduals,params0,args = (means_lsq, stds_lsq**2),full_output = 1)

    params0 = [0.2,40]

    k, kcov_x, info, mesg, success = leastsq(residuals,params0,args = (means_lsq, stds_lsq**2),full_output = 1)

    print 'pcov_x = {0}'.format(pcov_x)
    print 'kcov_x = {0}'.format(kcov_x)

    print 'Read noise = {0}'.format(np.sqrt(p[3]))
    print 'Gain = {0}'.format(1./k[0])

    x = np.arange(40000)
    y = np.arange(0,100)

    print 'Read noise'
    # Finding chi squared from the residuals and the errors
    chi_sqr = sum((polyresiduals(p,means_lsq,stds_lsq**2))**2 / (stds_lsq**2))
    print "Chi squared is", "%.4g" % chi_sqr

    dof = len(means_lsq)-len(p)

    # Finding the reduced chi squared
    red_chi_sqr = chi_sqr/(dof)
    print "Reduced chi squared is", "%.4g" % red_chi_sqr

    # And finding the residual variance
    res_var = np.sqrt(chi_sqr/(dof))
    print "Reduced variance is","%.4g" % res_var

    # Multiply to get parameter errors
    covs = pcov_x*res_var

    print 'Read noise error = {0}'.format(np.sqrt(covs[3][3]))

    print 'Gain'
    # Finding chi squared from the residuals and the errors
    chi_sqr = sum((residuals(p,means_lsq,stds_lsq**2))**2 / (stds_lsq**2))
    print "Chi squared is", "%.4g" % chi_sqr

    dof = len(means_lsq)-len(k)

    # Finding the reduced chi squared
    red_chi_sqr = chi_sqr/(dof)
    print "Reduced chi squared is", "%.4g" % red_chi_sqr

    # And finding the residual variance
    res_var = np.sqrt(chi_sqr/(dof))
    print "Reduced variance is","%.4g" % res_var

    # Multiply to get parameter errors
    covs = kcov_x*res_var

    print 'Slope error = {0}'.format(np.sqrt(covs[0][0]))
    print k

    plt.clf()
    plt.subplot(211)
    plt.plot(means,stds**2,'+')
    plt.xlabel('Pixel Mean')
    plt.ylabel('Pixel Variance')
    plt.title('Full 1024 Sample Set')
    plt.xlim(min(means),max(means))
    plt.subplot(212)
    plt.plot(x, guess(k,x),'r',linewidth = 2,label = 'Linear Fit')
    plt.plot(x, polyguess(p,x),'g:',linewidth = 5,label = 'Read Noise Fit')
    plt.plot(means_lsq, stds_lsq**2, 'bo',markersize = 1)
    plt.title('Mean and variation of each pixel in a selcted set of 200 samples')
    plt.xlabel('Pixel Mean')
    plt.ylabel('Pixel Variance')
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.show()

    plt.clf()
    plt.plot(np.log(means),np.log(stds),'o')
    plt.title('Photon Transfer Curve')
    plt.xlabel('Log of the mean signal')
    plt.ylabel('Log of the standard deviation of the signal')
    plt.show()

#f = import_files(files)
#d = import_files(dark_files)

#noise(f)

