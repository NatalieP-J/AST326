from PMT import *
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.misc import factorial
from scipy.stats import poisson

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
    #plt.plot(hr,hist,drawstyle = 'steps-mid')
    #plt.savefig('{0}.png'.format(somename))
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
    std = np.sqrt(sum(newlist)/len(hist)-1)
    #print 'The standard deviation is {0}'.format(std)
    #print 'The mean is {0}'.format(mean)
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
    for i in range(len(tsamp)):
        name = 'section7-' + time + '-{0}-{1}'.format(nsamp,tsamp[i])
        counts = np.loadtxt(name+'.dat')
        mean = sum(counts)/len(counts)
        means.append(mean)
        vals = []
        for j in range(len(counts)):
            val = ((counts[j])-mean)**2
            vals.append(val)
        std = np.sqrt(sum(vals)/(len(vals)-1))
        stds.append(std)
        var = std**2
        variance.append(var)
        total_counts.append(counts)
    name = 'section7-' + time +  'subplot.png'
    f = plt.figure()
    ax = plt.subplot(211)
    ax.set_title("Mean Count vs Variance on a Linear Scale") 
    ax.plot(means,variance,'o',label = 'mean vs variance')
    x = np.arange(min(means),max(means))
    y = x
    ax.plot(x,y,label = 'x=y')
    ax.set_xlabel('Mean Count')
    ax.set_ylabel('Variance')
    ax = plt.subplot(212)
    ax.set_title("Mean Count vs Variance on a Log Scale")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(means,variance,'o',label = 'mean vs variance')
    ax.plot(x,y, label = 'x=y')
    ax.set_xlabel('Mean Count on Log Scale')
    ax.set_ylabel('Variance on Log Scale')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(name)
    plt.show()
    return means,stds,variance,tsamp

#Section 8 - Poisson
def poisson2(x,mean):
    newlist = []
    for i in range(len(x)):
        factx = int(factorial(x[i]))
        p = (mean**x[i]/factx)*np.exp(-mean)
        newlist.append(p)
    return newlist
    
def gaussian(x,mean,std,maxval):
    newlist = []
    for k in range(len(x)):
        diff = x[k]-mean
        div = diff/std
        exp = np.exp(-0.5*(div)**2)
        const = 1./(std*np.sqrt(2*np.pi))
        g = const*exp
        newlist.append(g)
    return newlist

def section8new():
    tsamp = 0.05
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
names = ['section8-2013-9-23_17:1:54-0.003-1000',
            'section8-2013-9-23_17:1:13-0.003-1000',
            'section8-2013-9-23_16:14:8-0.003-1000',
            'section8-2013-9-23_17:0:2-0.003-1000',
            'section8-2013-9-23_16:9:28-0.003-1000',
            'section8-2013-9-23_16:18:55-0.003-1000']
            
names2 = ['section8-2013-9-23_18:23:59-0.3-1000',
    'section8-2013-9-23_18:30:20-0.3-1000',
    #'section8-2013-9-23_18:51:45-0.3-1000',
    'section8-2013-9-23_19:8:54-0.3-1000',
    'section8-2013-9-23_20:31:1-0.3-1000']
def section8(fnames):
    plt.figure()
    plt.suptitle('1000 0.3s Samples')
    for k in range(len(fnames)):
        counts = np.loadtxt(fnames[k]+'.dat')
        mean = np.mean(counts)
        sqrs = []
        for i in range(len(counts)):
            sqr = ((counts[i]-mean)**2)
            sqrs.append(sqr)
        std = np.sqrt((sum(sqrs))/len(sqrs)-1)
        name = fnames[k]
        h = histogram(counts,name)
        hr = h[0]
        step = 0.001
        x = np.arange(min(hr),max(hr),step)
        hist = h[1]
        prob = []
        total_hist = sum(hist)
        for i in range(len(hist)):
            p = float(hist[i])/total_hist
            prob.append(p)
        probmean = sum(prob)/len(prob)
        ax = plt.subplot(2,2,k)
        ax.plot(hr,prob,drawstyle = 'steps-mid',lw=1,label = 'Raw Data')
        ax.plot(hr,poisson(mean).pmf(hr),lw=1,label = 'Poisson Fit')
        ax.plot(x,gaussian(x,mean,std,max(prob)),lw=1,label = 'Gaussian Fit')
        ax.vlines(mean,0,max(prob),color='m',lw=2,label = 'Mean Value')
        #ax.set_ylim(0,0.30)
        #ax.vlines(mean-std,0,max(prob),color='r',lw=1)
        #ax.vlines(mean+std,0,max(prob),color='r',lw=1)
        #ax.vlines(mean+np.sqrt(mean),0,max(prob),color='g')
        #ax.vlines(mean-np.sqrt(mean),0,max(prob),color='g')
        ax.set_xlabel('Counts per Sample')
        ax.set_ylabel('Probability')
        g = gaussian(x,mean,std,max(prob))
        intg=[]
        print x
        for i in range(len(x)):
            rect = 0.01*g[i]
            intg.append(rect)
        print 'integration of gaussian:',sum(intg)
    plt.legend(loc='best',prop = {'size':8})
    plt.tight_layout()
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
            print 'working on {0} of {1}'.format(k,nsamps[j])
            counts = get_counts(tsamp,nsamps[j])
            save_data(name,counts)
            mean = float(sum(counts))/len(counts)

def section9(nmax,time):
    pts = []
    nsamps = []
    predict = []
    moms = []
    for i in range(1,nmax+1):
        n=2**i
        print n
        nsamps.append(n)
    tsamp = 0.01
    for j in range(len(nsamps)):
        means = []
        stds = []
        for k in range(15):
            sqrs = []
            name = 'section9-'+time+'-{0}-{1}_{2}'.format(tsamp,nsamps[j],k+1)
            data = np.loadtxt(name+'.dat')
            for d in range(len(data)):
                r = 0
                if data[d]==0:
                    np.delete(data,d)
                    r+=1
                if r>0:
                    print 'removed {0} outliers'.format(r)
            avg = float(sum(data))/len(data)
            means.append(avg)
            for i in range(len(data)):
                sqr = (float(data[i])-avg)**2
                sqrs.append(sqr)
            std = np.sqrt((sum(sqrs))/(len(sqrs)-1))
            stds.append(std)
        mom = float(sum(means))/len(means)
        moms.append(mom)
        totsqrs=[]
        val = float(sum(stds)/len(stds))
        print val
        for m in range(len(stds)):
            sqr = (stds[m]-val)**2
            totsqrs.append(sqr)
        print sum(totsqrs),len(totsqrs),j
        sdom = np.sqrt(sum(totsqrs)/(len(totsqrs)-1))
        sdom = val/np.sqrt(nsamps[j])
        point = [mom,sdom,nsamps[j]]
        pts.append(point)
        pre = []
        for i in range(len(means)):
            predict_std = np.sqrt(means[i])
            pre.append(predict_std)
        presqrs = []
        avg = sum(pre)/len(pre)
        predictsdom = avg/np.sqrt(nsamps[j])
        predict.append(predictsdom)
        print 'stds{0}='.format(j), stds
    print 'predict=',predict
    return pts

mean = [11.166666666666666,
    10.966666666666667,
    11.341666666666667,
    11.4125,
    11.53125,
    11.764583333333333,
    11.816666666666666,
    12.094010416666666,
    12.473152294035726,
    13.377417410977667]

std = [1.8986210534916925,
    1.0869714925491674,
    0.82198885932349863,
    0.66797303645313888,
    0.55878401974367553,
    0.32289562462724913,
    0.27113497308186602,
    0.2209633826330315,
    0.098486653654511619,
    0.11244532282960551]

nsamp = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

predict= [0.067488882675193762, 
        0.23781641715913096, 
        0.35340043300542923, 
        0.08099173787387802, 
        0.14547874551726847, 
        0.020731233035945144, 
        0.018453900220719983, 
        0.0046901457184241905, 
        0.052911022184930039, 
        0.11571141463783768]

predict = [2.3373915313460754, 1.6528907097169716, 1.1887096960872059, 0.8442058219535119, 0.60015477428651642, 0.42869671226694905, 0.30381201523892476, 0.21734787237850098, 0.15607661972824169, 0.11428625444222892]
std = [1.9666666666666663,
 1.5879349009259869,
 1.1409111661949083,
 0.86417719774040147,
 0.60247831977942778,
 0.43942516304591184,
 0.30582150408519387,
 0.22040346068516439,
 0.15823186984664031,
 0.11640214503365261]

f = plt.figure()
ax = plt.subplot(121)
ax.plot(nsamp,mean,'o',markersize=8,color='blue')
plt.ylabel('Mean of the Mean For 15 Trials')
plt.xlabel('Number of Samples') 
ax = plt.subplot(122)
ax.plot(nsamp,std,'o',color='orange',markersize=8,label='Standard Deviation \nof the Mean')
ax.plot(nsamp,predict,linewidth=2,color='green',markersize=8,label='Predicted Standard \nDeviation of the Mean')
ax.legend(loc='best',prop = {'size':12})
plt.ylabel('Standard Deviation of the Mean for 15 Trials')
plt.xlabel('Number of Samples') 
plt.tight_layout()
plt.show()
         
            
    




    
    
    
