
import numpy as np
import matplotlib.pyplot as plt

def avg_per_cycle(p,metric):
    # Read the data file associated with the given metric for the run with p processes
    # return the average of the metric over all cycles
    fm = open ('P'+str(p)+'-' + metric + '.data','r')
    cycle=0
    for line in fm:
        cycle = cycle + 1
        value = float(line)
    return value / cycle

def max_per_cycle(p,metric):
    # Read the data file associated with the given metric for the run with p processes
    # return the maximum of the metric over all cycles
    fm = open ('P'+str(p)+'-' + metric + '.data','r')
    max_value = 0
    for line in fm:
        value = float(line)
        if (value > max_value):
            max_value = value
    return max_value

pf = open ('procs.data','r')

procs =          []
bytes_highest =  []
time_adapt =     []
time_compute =   []
time_cycle =     []
time_initial =   []
time_output =    []
time_prepare =   []
time_refresh =   []
time_wallclock = []

x=[]
xp=[]
count = 0

for line in pf:

    p = int(line);
    procs.append(line)
    x.append(count)
    xp.append(count+0.4)
    count = count + 1

    bytes_highest.append  ( max_per_cycle(p,'bytes-highest') )
    time_adapt.append     ( avg_per_cycle(p,'time-adapt') )
    time_compute.append   ( avg_per_cycle(p,'time-compute') )
    time_cycle.append     ( avg_per_cycle(p,'time-cycle') )
    time_initial.append   ( avg_per_cycle(p,'time-initial') )
    time_output.append    ( avg_per_cycle(p,'time-output') )
    time_prepare.append   ( avg_per_cycle(p,'time-prepare') )
    time_refresh.append   ( avg_per_cycle(p,'time-refresh') )
    time_wallclock.append ( avg_per_cycle(p,'wallclock') )

pf.close()

print p,' time_cycle  = ',time_cycle

print p,' time_compute  = ',time_compute
print p,' time_adapt  = ',time_adapt
print p,' time_refresh  = ',time_refresh
print p,' time_prepare  = ',time_prepare

t = ()
print t

plot_wallclock = plt.bar(x,time_wallclock,color='r')

plot_cycle     = plt.bar(x,time_cycle,color='g')

plot_compute   = plt.bar(x,time_compute,color='b')
top = time_compute
plot_adapt     = plt.bar(x,time_adapt,  color='c',bottom=top)
top = [sum(pair) for pair in zip(top,time_adapt)]
plot_refresh   = plt.bar(x,time_refresh,color='m',bottom=top)
top = [sum(pair) for pair in zip(top,time_refresh)]
plot_prepare   = plt.bar(x,time_prepare,color='y',bottom=top)

plt.title("Enzo-P Scaling")
plt.legend( (plot_wallclock[0], 
             plot_cycle[0],
             plot_compute[0],
             plot_adapt[0],
             plot_refresh[0],
             plot_prepare[0]),
            ( 'wallclock', 'cycle','compute','adapt','refresh','prepare'),
            loc="upper left")
plt.xticks(xp,procs)
plt.grid(True)

plt.show()


#P16-bytes-highest.data	P16-time-cycle.data    P16-time-prepare.data
#P16-time-adapt.data	P16-time-initial.data  P16-time-refresh.data
#P16-time-compute.data	P16-time-output.data   P16-wallclock.data

# Data: 
#   Wallclock
#   cycle
#   compute
#   adapt
#   refresh
#   prepare
#   initial
#   output
