#!/usr/bin/python

from numpy import *
import numpy as np
import matplotlib.pyplot as plt

color = ['r','g','b','y','c','m','k']

file = [
    'adapt',
    'compute',
    'cycle',
    'initial',
    'output',
    'refresh',
    'simulation',
    'stopping']

plt.title("Time per Cycle")
plt.xlabel("Cycle")
plt.ylabel("Time (s)")

# plt.text(10,5.0e9,"GNU")
# plt.text(10,4.5e9,"Intel")
# plt.text(10,4.0e9,"PGI")

# plt.plot(8,5.08e9,'ro',label="GNU")
# plt.plot(8,4.58e9,'go',label="Intel")
# plt.plot(8,4.08e9,'bo',label="PGI")

for i in range(len(file)):

   f = loadtxt('usec-'+file[i]+'.data',dtype=int)
   plt.plot(f[:,0],10e-6*f[:,1], label=file[i])

plt.legend(loc="upper left")

plt.grid(True)

print ("Writing 'usec.pdf'")

plt.savefig("usec.pdf", format='pdf')
