#!/usr/bin/python

from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import glob
import os.path
import string

def open_plot(plt,title,label_x,label_y):

    plt.clf()

    plt.axes([0.1,0.1,0.8,0.8])
    plt.figtext(0.5,0.95,title,fontsize=18, ha='center')
    # plt.figtext(0.5,0.915, "N7 P(1,*,1) net-gcc", fontsize=12, ha='center')
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.grid(True)

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

open_plot(plt,"Enzo-E: Cumulative time in adapt","cycle","time (s)");

f_adapt =  loadtxt('adapt.data',dtype=float)
adapt_x_total = f_adapt[:,0]
adapt_y_total = f_adapt[:,1]
plt.plot(adapt_x_total,1e-6*adapt_y_total, label='adapt', color=c[0],marker=m[0], ls=l, lw=w)

i=1
for file in list(glob.glob('adapt_*data')):
    if (not file.endswith(("_post.data"))):
        f_adapt =  loadtxt(file,dtype=float)
        x_cycle = f_adapt[:,0]
        y_adapt = f_adapt[:,1]
        name=os.path.splitext(os.path.basename(file))[0]
        plt.plot(x_cycle,1e-6*y_adapt, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
        i=i+1

plt.legend(loc='lower center',ncols=4)
plt.yscale("log") 
print ("Writing 'plot_adapt'")
plt.savefig("plot_adapt.pdf", format='pdf')
plt.savefig("plot_adapt.png", format='png')

# ======================================================================

open_plot(plt,"Enzo-E: Cumulative time in refresh","cycle","time (s)");

f_refresh =  loadtxt('refresh.data',dtype=float)
refresh_x_total = f_refresh[:,0]
refresh_y_total = f_refresh[:,1]
plt.plot(refresh_x_total,1e-6*refresh_y_total, label='refresh', color=c[0],marker=m[0], ls=l, lw=w)

i=1
for file in list(glob.glob('refresh_*data')):
    if (not file.endswith(("_post.data"))):
        f_refresh =  loadtxt(file,dtype=float)
        x_cycle = f_refresh[:,0]
        y_refresh = f_refresh[:,1]
        name=os.path.splitext(os.path.basename(file))[0]
        plt.plot(x_cycle,1e-6*y_refresh, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
        i=i+1

plt.legend(loc='lower center',ncols=4)
plt.yscale("log") 
print ("Writing 'plot_refresh'")
plt.savefig("plot_refresh.pdf", format='pdf')
plt.savefig("plot_refresh.png", format='png')

# ======================================================================

open_plot(plt,"Enzo-E: Cumulative time in method","cycle","time (s)");

f_method =  loadtxt('method.data',dtype=float)
method_x_total = f_method[:,0]
method_y_total = f_method[:,1]
plt.plot(method_x_total,1e-6*method_y_total, label='method', color=c[0],marker=m[0], ls=l, lw=w)

i=1
for file in list(glob.glob('method_*data')):
    if (not file.endswith(("_post.data"))):
        f_method =  loadtxt(file,dtype=float)
        x_cycle = f_method[:,0]
        y_method = f_method[:,1]
        name=os.path.splitext(os.path.basename(file))[0]
        plt.plot(x_cycle,1e-6*y_method, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
        i=i+1

plt.legend(loc='lower center',ncols=4)
plt.yscale("log") 
print ("Writing 'plot_method'")
plt.savefig("plot_method.pdf", format='pdf')
plt.savefig("plot_method.png", format='png')

# ======================================================================

open_plot(plt,"Enzo-E: Cumulative time in solver","cycle","time (s)");

f_solver =  loadtxt('solver.data',dtype=float)
solver_x_total = f_solver[:,0]
solver_y_total = f_solver[:,1]
plt.plot(solver_x_total,1e-6*solver_y_total, label='solver', color=c[0],marker=m[0], ls=l, lw=w)

i=1
for file in list(glob.glob('solver_*data')):
    if (not file.endswith(("_post.data"))):
        f_solver =  loadtxt(file,dtype=float)
        x_cycle = f_solver[:,0]
        y_solver = f_solver[:,1]
        name=os.path.splitext(os.path.basename(file))[0]
        plt.plot(x_cycle,1e-6*y_solver, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
        i=i+1

plt.legend(loc='lower center',ncols=4)
plt.yscale("log") 
print ("Writing 'plot_solver'")
plt.savefig("plot_solver.pdf", format='pdf')
plt.savefig("plot_solver.png", format='png')

# ======================================================================

open_plot(plt,"Enzo-E: Cumulative times","cycle","time (s)");

f_cycle =  loadtxt('cycle.data',dtype=float)
cycle_x_total = f_cycle[:,0]
cycle_y_total = f_cycle[:,1]
i=0
plt.plot(cycle_x_total,cycle_y_total, label='cycle', color=c[i],marker=m[i], ls=l, lw=w)
i=i+1
plt.plot(method_x_total,1e-6*method_y_total, label='method', color=c[i],marker=m[i], ls=l, lw=w)
i=i+1
plt.plot(solver_x_total,1e-6*solver_y_total, label='solver', color=c[i],marker=m[i], ls=l, lw=w)
i=i+1
plt.plot(refresh_x_total,1e-6*refresh_y_total, label='refresh', color=c[i],marker=m[i], ls=l, lw=w)
i=i+1
plt.plot(adapt_x_total,1e-6*adapt_y_total, label='adapt', color=c[i],marker=m[i], ls=l, lw=w)
i=i+1

plt.legend(loc='lower center',ncols=4)
plt.yscale("log") 
print ("Writing 'plot_total'")
plt.savefig("plot_total.pdf", format='pdf')
plt.savefig("plot_total.png", format='png')
