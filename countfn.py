from PMT import *
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def get_counts(tsamp,nsamp):
    counts = photoncount(tsamp,nsamp)
    return counts

def mean_rate(somelist,int_time):
    newlist=[]
    length = len(somelist)
    for i in range(length):
        rate = somelist[i]/int_time
        newlist.append(rate)
    avg_rate = sum(newlist)/len(newlist)
    return avg_rate

def save_data(data_name,somelist):
    np.savetxt('{0}.dat'.format(data_name),somelist,fmt='%i')

def standard_deviation(somearray):
    dat = [] 
    avg = np.sum(somearray)/(np.size(somearray))
    sqravg = np.sum(somearray*somearray)/np.size(x)
    std = np.sqrt(sqravg - avg*avg)
    print 'The mean is {0}'.format(avg)
    print 'The standard deviation is {0}'.format(std)
    dat.append(avg)
    dat.append(std)
    return dat

def histogram(somearray,somename):
    hmin = min(somearray)
    hmax = max(somearray)
    hr = np.arange(hmin,hmax+1)
    hist = np.array([np.where(somearray == i)[0].size for i in hr])
    plt.plot(hr,hist,drawstyle = 'steps-mid')
    plt.show()
    plt.savefig('{0}.png'.format(somename))

def pmt_loop (tsamp,nsamp,max_i):
    raw_data = []
    for i in range(max_i):
        name = '-{0}_{1}_{2}'.format(tsamp,nsamp,max_i)
        stamp = datetime.now()
        year = stamp.year
        month = stamp.month
        day = stamp.day
        hour = stamp.hour
        minute = stamp.minute
        second = stamp.second
        timestr = '{0}-{1}-{2}_{3}:{4}:{5}'.format(year,month,day,hour,minute,second)
        name = timestr + name
        count = get_counts(tsamp,nsamp)
        save_data(name,count)
        histogram(count,name)
        raw_data.append(count)
    return raw_data
