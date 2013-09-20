from PMT import *
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.misc import factorial

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

def histogram(somearray,somename):
    hmin = min(somearray)
    hmax = max(somearray)
    hr = np.arange(hmin,hmax+1)
    hist = np.array([np.where(somearray == i)[0].size for i in hr])
    plt.plot(hr,hist,drawstyle = 'steps-mid')
    plt.savefig('{0}.png'.format(somename))
    return hr,hist

def plot_raw(somearray,units):
    plt.plot(somearray,drawstyle = 'steps-mid')
    plt.xlabel('Time ({0})'.format(units))
    plt.ylabel('Number of Counts')
    plt.show()

def timestamp():
    stamp = datetime.now()
    year = stamp.year
    month = stamp.month
    day = stamp.day
    hour = stamp.hour
    minute = stamp.minute
    second = stamp.second
    timestr = '{0}-{1}-{2}_{3}:{4}:{5}'.format(year,month,day,hour,minute,second)
    return timestr

def standard_deviation(somearray):
    dat = []
    max_val = max(somearray)
    min_val = min(somearray)
    vals = np.arange(min_val,max_val+1)
    hist = np.array([np.where(somearray == i)[0].size for i in vals])
    mean = sum(hist)/len(hist)
    for i in range(len(hist)):
        newlist = []
        sqr = (hist[i] - mean)**2
        newlist.append(sqr)
    std = np.sqrt(sum(newlist)/len(hist))
    print 'The standard deviation is {0}'.format(std)
    print 'The mean is {0}'.format(mean)
    dat.append(mean)
    dat.append(std)
    return dat

def pmt_loop (tsamp,nsamp,max_i):
    raw_data = []
    timestr = timestamp()
    for i in range(max_i):
        print 'Working on trial number {0}'.format(i)
        name = '-{0}_{1}_{2}'.format(tsamp,nsamp,i)
        name = timestr + name
        count = get_counts(tsamp,nsamp)
        save_data(name,count)
        histogram(count,name)
        raw_data.append(count)
    for i in range(len(raw_data)):
        if tsamp == 0.001:
            units = 'microseconds'
        if tsamp == 0.01:
            units = 'milliseconds'
        if tsamp == 1:
            units = 'seconds'
        #plot_raw(raw_data[i],units) 
    return raw_data


#Section 7
def section7new():
    total_counts=[]
    means = []
    stds = []
    variance = []
    tsamp = [0.001,0.005,0.01,0.05,0.1,0.5]
    nsamp = 100
    time = timestamp()
    for i in range(len(tsamp)):
        name = 'section7-' + time + '-{0}-{1}'.format(nsamp,tsamp[i])
        counts = get_counts(tsamp[i],nsamp)
        print counts
        save_data(name,counts)
        mean = sum(counts)/len(counts)
        means.append(mean)
        vals = []
        for j in range(len(counts)):
            val = (counts[j]-mean)**2
            vals.append(val)
        std = np.sqrt(sum(vals)/(len(vals)-1))
        stds.append(std)
        var = std**2
        variance.append(var)
        total_counts.append(counts)
    name = 'section7-' + time + '.png'
    f = plt.figure()
    plt.plot(means,variance,'-o')
    x = np.arange(min(means),max(means))
    y = x
    plt.plot(x,y)
    plt.savefig(name)
    plt.show()
    f = plt.figure()
    name = 'section7-' + time + '-logplot.png'
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(means,variance,'-o')
    plt.plot(x,y,'-o')
    plt.savefig(name)
    plt.show()
    return means,stds,variance,tsamp

def section7(time):
    total_counts=[]
    means = []
    stds = []
    variance = []
    tsamp = [0.001,0.005,0.01,0.05,0.1,0.5]
    nsamp = 100
    time = timestamp()
    for i in range(len(tsamp)):
        name = 'section7-' + time + '-{0}-{1}'.format(nsamp,tsamp[i])
        counts = np.loadtxt(name+'.dat')
        mean = sum(counts)/len(counts)
        means.append(mean)
        vals = []
        for j in range(len(counts)):
            val = (counts[j]-mean)**2
            vals.append(val)
        std = np.sqrt(sum(vals)/(len(vals)-1))
        stds.append(std)
        var = std**2
        variance.append(var)
        total_counts.append(counts)
    name = 'section7-' + time + '.png'
    f = plt.figure()
    plt.plot(means,variance,'-o')
    x = np.arange(min(means),max(means))
    y = x
    plt.plot(x,y)
    plt.savefig(name)
    plt.show()
    f = plt.figure()
    name = 'section7-' + time + '-logplot.png'
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(means,variance,'-o')
    plt.plot(x,y,'-o')
    plt.savefig(name)
    plt.show()
    return means,stds,variance,tsamp

#Section 8 - Poisson
def poisson(x,mean):
    newlist = []
    for i in range(len(x)):
        factx = int(factorial(x[i]))
        p = (mean**x[i]/factx)*np.exp(-mean)
        newlist.append(p)
    return newlist

def section8new():
    tsamp = 0.003
    nsamp = 1000
    time = timestamp()
    name = 'section8-' + time + '-{0}-{1}'.format(tsamp,nsamp)
    counts = get_counts(tsamp,nsamp)
    save_data(name,counts)
    mean = sum(counts)/len(counts)
    h = histogram(counts,name)
    hr = h[0]
    hist = h[1]
    prob = []
    total_hist = sum(hist)
    for i in range(len(hist)):
        p = float(hist[i])/total_hist
        prob.append(p)
    plt.figure()
    plt.plot(hr,prob,'-o')
    plt.plot(hr,poisson(hr,mean),'-o')
    plt.show()

def section8(fname):
    counts = np.loadtxt(fname+'.dat')
    mean = sum(counts)/len(counts)
    h = histogram(counts,name)
    hr = h[0]
    hist = h[1]
    prob = []
    total_hist = sum(hist)
    for i in range(len(hist)):
        p = float(hist[i])/total_hist
        prob.append(p)
    plt.figure()
    plt.plot(hr,prob,'-o')
    plt.plot(hr,poisson(hr,mean),'-o')
    plt.show()
    
#Section 9
def section9new(nmax):
    moms = []
    sdoms = []
    means = []
    stds = []
    nsamps = []
    time = timestamp()
    for i in range(1,nmax+1):
        n = 2**i
        nsamps.append(n)
    tsamp = 0.01
    for j in range(len(nsamps)):
        mean = []
        std = []
        for k in range(15):
            name = 'section9-' + time + '-{0}-{1}_{2}'.format(tsamp,nsamps[j],k+1)
            counts = get_counts(tsamp,nsamps[j])
            save_data(name,counts)
            mean = float(sum(counts))/len(counts)
    return counts



    
    
    
