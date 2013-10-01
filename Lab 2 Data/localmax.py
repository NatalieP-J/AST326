import matplotlib.pyplot as plt
import numpy as np

files = /h/ungrad1/stostad/lab2/neon/neon000

def avg(files):
    for indx in range(48):
        file1 = files + str(round(indx,2)) + 

def maxima(file, threshold = 2000):
    data = np.genfromtxt(file, skip_header = 17, skip_footer = 1)
    for i in range(len(data[1:-1])):
        if data[i,1] > data[i-1,1] and data[i,1] > data[i+1,1] and data[i,1] > threshold:
            print data[i]
