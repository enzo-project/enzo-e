#!/usr/bin/python

from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import glob
import os.path
import string

#--------------------------------------------------
def plot_open(plt,title,label_x,label_y):

    plt.clf()

    plt.axes([0.1,0.1,0.8,0.8])
    plt.figtext(0.5,0.95,title,fontsize=18, ha='center')
    # plt.figtext(0.5,0.915, "N7 P(1,*,1) net-gcc", fontsize=12, ha='center')
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.grid(True)

#--------------------------------------------------
def plot_total(plt,name,label,i=0,scale=1e-6):
    data =  loadtxt(name,dtype=float)
    total_x = data[:,0]
    total_y = data[:,1]
    plt.plot(total_x,scale*total_y, label=label, color=c[i%7],marker=m[i%7], ls=l, lw=w)
    return [total_x, total_y]

#--------------------------------------------------
def plot_list(plt,file_list,scale=1e-6):
    file_times={}
    # sort file names by revert end time
    for file in file_list:
        if (not file.endswith(("_post.data"))):
            data =  loadtxt(file,dtype=float)
            glob_x = data[:,0]
            glob_y = data[:,1]
            file_times[glob_y[len(glob_y)-1]] = file
    i=1
    # plot by revert end time
    for size in sorted(file_times,reverse=True):
        file=file_times[size]
        data =  loadtxt(file,dtype=float)
        glob_x = data[:,0]
        glob_y = data[:,1]
        name=os.path.splitext(os.path.basename(file))[0]
        plt.plot(glob_x,scale*glob_y, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
        i=i+1

def plot_write(name):
    print ("Writing 'plot_%s'" % (name))
    plt.savefig(("plot_%s.pdf") % (name), format='pdf')
    plt.savefig(("plot_%s.png") % (name), format='png')

# ======================================================================

c=['b', 'g', 'r', 'c', 'm', 'y', 'k']
m=['>', '^', '<', 'v', 'o', '+', 'x', 's' ]
m=['none', 'none', 'none', 'none', 'none', 'none', 'none', 'none']
l='-'
w=1

figure(figsize=(8,6), dpi=100)

#  adapt*data

#  refresh*data
#  method*data
#  solver*data
#  memory*data
#    [bytes*data]
#  mesh*data
#    [*blocks*data]
#  summary
#    cycle
#    method
#    solver
#    adapt
#    refresh

# ======================================================================
# ADAPT
# ======================================================================

plot_open(plt,"Enzo-E: Cumulative time in adapt","cycle","time (s)");
[adapt_x_total, adapt_y_total] = plot_total(plt,'adapt.data','adapt')
plot_list(plt,glob.glob('adapt_*data'))

plt.legend(loc='lower center',ncols=4)
plt.yscale("log")
plot_write('plot_adapt')

# ======================================================================

plot_open(plt,"Enzo-E: Cumulative time in refresh","cycle","time (s)");
[refresh_x_total, refresh_y_total] = plot_total(plt,'refresh.data','refresh')
plot_list(plt,glob.glob('refresh_*data'))

plt.legend(loc='lower center',ncols=4)
plt.yscale("log")
plot_write('plot_refresh')

# ======================================================================

plot_open(plt,"Enzo-E: Cumulative time in method","cycle","time (s)");
[method_x_total, method_y_total] = plot_total(plt,'method.data','method')
plot_list(plt,glob.glob('method_*data'))

plt.legend(loc='lower center',ncols=4)
plt.yscale("log")
plot_write('plot_method')

# ======================================================================

plot_open(plt,"Enzo-E: Cumulative time in solver","cycle","time (s)");
[solver_x_total, solver_y_total] = plot_total(plt,'solver.data','solver')
plot_list(plt,glob.glob('solver_*data'))

plt.legend(loc='lower center',ncols=4)
plt.yscale("log")
plot_write('plot_solver')

# ======================================================================

plot_open(plt,"Enzo-E: Cumulative times","cycle","time (s)");

plot_total(plt,'cycle.data','cycle',0,1.0)
plot_list(plt,["method.data", "solver.data", "refresh.data", "adapt.data"])

plt.legend(loc='lower center',ncols=4)
plt.yscale("log")
plot_write('plot_total')
