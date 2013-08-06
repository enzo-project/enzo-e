
import numpy as np
import matplotlib.pyplot as plt


dir = '.'

files = [dir + '/time-startup.data',
         dir + '/time-total.data',
         dir + '/procs.data']

data = [[],[],[]]

for i in range(len(files)):
   f = open(files[i],'r')
   for line in f:
      data[i] = data[i] + [float(line)]
   f.close()

plt.title("(DRAFT) Sedov Blast Array Weak Scaling\nstartup time and time per cycle")
plt.xlabel("Processors $P$")
plt.ylabel("Time T_P (s)")

time_startup = np.array(data[0])
time_total   = np.array(data[1])
procs        = np.array(data[2])

x = procs
y1 = time_startup
y2 = time_total / 10

plt.plot(x,y1,'ro-', label="startup time")
plt.plot(x,y2,'go-', label="time per cycle")
# a = plt.gca()
# a.set_ylim([0,1.1])
plt.xscale('log')
plt.yscale('log')

plt.legend(loc="upper left")

plt.grid(True)

plt.savefig("time-total.pdf", format='pdf')

