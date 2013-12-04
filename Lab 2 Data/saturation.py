import matplotlib.pyplot as plt
import numpy as np


files = 'saturation2/sat2_00'

def sat(files, nofiles = 501):
    all_data = np.zeros((nofiles, 2048, 2))
    # First loop over the files to get data from each and put the data
    # in all_data:
    for indx in range(nofiles):
        file1 = files + str(indx).zfill(len(str(nofiles))) + '.txt'
        data1 = np.genfromtxt(file1, skip_header = 17, skip_footer = 1)
        all_data[indx] = data1
    datamean = np.zeros((2,2048))
    plt.clf()
    for i in [100, 500, 1000, 1250, 1500, 1750, 2000]:
        plt.plot(all_data[:,i,1]/(np.median(all_data[:,i,1])))
    # plt.subplot(321)
    # plt.plot(all_data[:,100,1])
    # plt.subplot(322)
    # plt.plot(all_data[:,600,1])
    # plt.subplot(323)
    # plt.plot(all_data[:,1000,1])
    # plt.subplot(324)
    # plt.plot(all_data[:,1300,1])
    # plt.subplot(325)
    # plt.plot(all_data[:,1600,1])
    # plt.subplot(326)
    # plt.plot(all_data[:,2000,1])
    plt.show()
    return all_data

all = sat(files)
a = [100, 500, 1000, 1250, 1500, 1750, 2000]
for i in a:
    print all[260, i, 1]
