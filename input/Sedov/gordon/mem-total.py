
import numpy as np
import matplotlib.pyplot as plt


dir = '.'

files = [dir + '/bytes-highest.data',
         dir + '/procs.data']

data = [[],[]]

for i in range(len(files)):
   f = open(files[i],'r')
   for line in f:
      data[i] = data[i] + [float(line)]
   f.close()

plt.title("(DRAFT) Sedov Blast Array Weak Scaling\naverage heap use per process")
plt.xlabel("Processors $P$")
plt.ylabel("MBytes / process")

bytes_highest = np.array(data[0])
procs         = np.array(data[1])

x = procs
y1 = bytes_highest*1e-6

plt.plot(x,y1,'ro-', label="average heap use per process")
plt.ylim([0,500])

plt.xscale('log')

plt.legend(loc="lower right")

plt.grid(True)


plt.savefig("mem-total.pdf", format='pdf')

