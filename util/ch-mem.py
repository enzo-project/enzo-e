#!/usr/bin/python

from numpy import *
import numpy as np
import matplotlib.pyplot as plt

color = ['r','g','b','y','c','m']

file = [
    'adapt',
    'refresh',
    'compute',
    'output',
    'stopping']

plt.title("Bytes per Block")
plt.xlabel("Block")
plt.ylabel("Bytes")

# plt.text(10,5.0e9,"GNU")
# plt.text(10,4.5e9,"Intel")
# plt.text(10,4.0e9,"PGI")

# plt.plot(8,5.08e9,'ro',label="GNU")
# plt.plot(8,4.58e9,'go',label="Intel")
# plt.plot(8,4.08e9,'bo',label="PGI")

for i in range(len(file)):

   f = loadtxt('mem-'+file[i]+'.data',dtype=int)
   plt.semilogy(f[:,0],1e-6*f[:,1], label=file[i])

plt.legend(loc="upper right")

plt.grid(True)

print ("Writing 'mem.pdf'")

plt.savefig("mem.pdf", format='pdf')
