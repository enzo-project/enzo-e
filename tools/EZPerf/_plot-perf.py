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
    # plt.figtext(0.5,0.915, 'N7 P(1,*,1) net-gcc', fontsize=12, ha='center')
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.grid(True)

#--------------------------------------------------
def plot_total(plt,filename,label,type='total',scale=1e-6):
    data =  loadtxt(filename,dtype=float)
    total_x = data[:,0]
    total_y = data[:,1]
    if (type == 'cycle'):
        n=len(total_x)
        total_x = total_x[1:n-1]
        total_y = total_y[2:n] - total_y[1:n-1]
    plt.plot(total_x,scale*total_y, label=label, color=c[0],marker=m[0], ls=l, lw=w)
    return [total_x, total_y]

#--------------------------------------------------
def plot_list(plt,file_list,type='total',scale=1e-6,sort=True):
    if (sort):
        # sort file names by revert end time
        file_times = []
        for file in file_list:
            if (not file.endswith(('_post.data'))):
                data =  loadtxt(file,dtype=float)
                glob_x = data[:,0]
                glob_y = data[:,1]
                print (len(glob_y)-1,file)
                file_times.append([glob_y[len(glob_y)-1],file])
        i=1
        # plot by revert end time
        for row in sorted(file_times,reverse=True):
            file = row[1]
            data =  loadtxt(file,dtype=float)
            glob_x = data[:,0]
            glob_y = data[:,1]
            if (type == 'cycle'):
                n=len(glob_x)
                glob_x = glob_x[1:n-1]
                glob_y = glob_y[2:n] - glob_y[1:n-1]
            name=os.path.splitext(os.path.basename(file))[0]
            plt.plot(glob_x,scale*glob_y, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
            i=i+1
    else:
        i=1
        # plot by revert end time
        for file in sorted(file_list):
            print ("plot_list unsorted ",file)
            data =  loadtxt(file,dtype=float)
            glob_x = data[:,0]
            glob_y = data[:,1]
            if (type == 'cycle'):
                n=len(glob_x)
                glob_x = glob_x[1:n-1]
                glob_y = glob_y[2:n] - glob_y[1:n-1]
            name=os.path.splitext(os.path.basename(file))[0]
            plt.plot(glob_x,scale*glob_y, label=name, color=c[i%7],marker=m[i%7], ls=l, lw=w)
            i=i+1

#--------------------------------------------------

def plot_write(name,html):
    print ("Writing '%s'" % (name))
    plt.savefig(('%s.pdf') % (name), format='pdf')
    plt.savefig(('%s.png') % (name), format='png')
    html_image(html,name + ".png")

#====================================================================== 

def html_start():

    html=open('index.html', 'w')
    html.write('<HTML>\n')
    html.write('  <HEAD>\n')
    html.write('    <link href="cello.css" rel="stylesheet" type="text/css">\n')
    html.write('  </head>\n')
    html.write('    <title>Enzo-E / Cello EZPerf</title>\n')
    html.write('    <body>\n')
    html.write('      <h1>Enzo-E / Cello EZPerf</h1>\n')
    return html

def html_stop(html):
    html.write('    </body>\n')
    html.write('</html>\n')

def html_start_table(html):
    html.write('      <table>\n')

def html_stop_table(html):
    html.write('      </table>\n')

def html_start_row(html):
    html.write('         <tr>\n')

def html_stop_row(html):
    html.write('         </tr>\n')

def html_image(html,image):
    html.write('            <td>\n')
    html.write('              <a href="' +image+'"><img width=512 src="'+image+'"></img></a>\n')
    html.write('            </td>\n')

# ----------------------------------------------------------------------

c=['b', 'g', 'r', 'c', 'm', 'y', 'k']
m=['>', '^', '<', 'v', 'o', '+', 'x', 's' ]
m=['none', 'none', 'none', 'none', 'none', 'none', 'none', 'none']
l='-'
w=1

figure(figsize=(8,6), dpi=100)


html=html_start()

html_start_table(html)

# ======================================================================
# ROW 1: TIMES SUMMARY
# ======================================================================


html_start_row(html)

# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: cumulative times','cycle','time (s)');
plot_total(plt,'cycle.data','cycle',scale=1.0)
if os.path.exists('smp.data'):
    plot_list(plt,['method.data', 'solver.data', 'refresh.data', 'adapt.data', 'smp.data'])
else:
    plot_list(plt,['method.data', 'solver.data', 'refresh.data', 'adapt.data'])
plt.legend(loc='lower center',ncols=3)
plt.yscale('log')
plot_write('plot_time_total',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: per-cycle times','cycle','time (s)');
plot_total(plt,'cycle.data','cycle',scale=1.0,type='cycle')
if os.path.exists('smp.data'):
    plot_list(plt,['method.data', 'solver.data', 'refresh.data', 'adapt.data', 'smp.data'],type='cycle')
else:
    plot_list(plt,['method.data', 'solver.data', 'refresh.data', 'adapt.data'],type='cycle')
plt.legend(loc='lower center',ncols=3)
plot_write('plot_time_cycle',html)
# ----------------------------------------------------------------------

html_stop_row(html)

# ======================================================================
# ROW 2: MEMORY, BLOCKS, BALANCE
# ======================================================================

html_start_row(html)

# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: memory usage','cycle','bytes');
plot_list(plt,glob.glob('memory_*data'))
plt.legend(loc='lower center',ncols=3)
ym,yp = plt.ylim()
plt.ylim(0,yp)
plot_write('plot_memory',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: blocks per level','cycle','blocks');
plot_total(plt,'mesh_total-blocks.data','mesh_blocks-total',scale=1.0)
plot_list(plt,glob.glob('mesh_blocks*.data'),scale=1.0,sort=False)
plt.legend(loc='upper left',ncols=1)
plot_write('plot_mesh',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: load-balance efficiency','cycle','efficiency');
plot_list(plt,glob.glob('balance_eff-*data'),scale=1.0)
plt.legend(loc='lower center',ncols=3)
plt.ylim(0,1)
plot_write('plot_balance_eff',html)
# ----------------------------------------------------------------------


html_stop_row(html)

# ======================================================================
# ROW 3: TOTAL TIMES
# ======================================================================

html_start_row(html)

# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: cumulative time in method','cycle','time (s)');
[method_x_total, method_y_total] = plot_total(plt,'method.data','method')
plot_list(plt,glob.glob('method_*data'))
plt.legend(loc='lower center',ncols=3)
plt.yscale('log')
plot_write('plot_method_total',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: cumulative time in solver','cycle','time (s)');
[solver_x_total, solver_y_total] = plot_total(plt,'solver.data','solver')
plot_list(plt,glob.glob('solver_*data'))
plt.legend(loc='lower center',ncols=3)
plt.yscale('log')
plot_write('plot_solver_total',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: cumulative time in adapt','cycle','time (s)');
[adapt_x_total, adapt_y_total] = plot_total(plt,'adapt.data','adapt')
plot_list(plt,glob.glob('adapt_*data'))
plt.legend(loc='lower center',ncols=3)
plt.yscale('log')
plot_write('plot_adapt_total',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: cumulative time in refresh','cycle','time (s)');
[refresh_x_total, refresh_y_total] = plot_total(plt,'refresh.data','refresh')
plot_list(plt,glob.glob('refresh_*data'))
plt.legend(loc='lower center',ncols=3)
plt.yscale('log')
plot_write('plot_refresh_total',html)
# ----------------------------------------------------------------------
if os.path.exists('smp.data'):
    plot_open(plt,'Enzo-E: cumulative time in SMP','cycle','time (s)');
    [smp_x_total, smp_y_total] = plot_total(plt,'smp.data','smp')
    plot_list(plt,glob.glob('smp_*data'))
    plt.legend(loc='lower center',ncols=3)
    plt.yscale('log')
    plot_write('plot_smp_total',html)
# ----------------------------------------------------------------------

html_stop_row(html)

# ======================================================================
# ROW 4: CYCLE TIMES
# ======================================================================

html_start_row(html)

# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: per-cycle time in method','cycle','time (s)');
[method_x_total, method_y_total] = plot_total(plt,'method.data','method',type='cycle')
plot_list(plt,glob.glob('method_*data'),type='cycle')
plt.legend(loc='lower center',ncols=3)
plot_write('plot_method_cycle',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: per-cycle time in solver','cycle','time (s)');
[solver_x_total, solver_y_total] = plot_total(plt,'solver.data','solver',type='cycle')
plot_list(plt,glob.glob('solver_*data'),type='cycle')
plt.legend(loc='lower center',ncols=3)
plot_write('plot_solver_cycle',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: per-cycle time in adapt','cycle','time (s)');
[adapt_x_total, adapt_y_total] = plot_total(plt,'adapt.data','adapt',type='cycle')
plot_list(plt,glob.glob('adapt_*data'),type='cycle')
plt.legend(loc='lower center',ncols=3)
plot_write('plot_adapt_cycle',html)
# ----------------------------------------------------------------------
plot_open(plt,'Enzo-E: per-cycle time in refresh','cycle','time (s)');
[refresh_x_total, refresh_y_total] = plot_total(plt,'refresh.data','refresh',type='cycle')
plot_list(plt,glob.glob('refresh_*data'),type='cycle')
plt.legend(loc='lower center',ncols=3)
plot_write('plot_refresh_cycle',html)
# ----------------------------------------------------------------------
if os.path.exists('smp.data'):
    plot_open(plt,'Enzo-E: per-cycle time in smp','cycle','time (s)');
    [smp_x_total, smp_y_total] = plot_total(plt,'smp.data','smp',type='cycle')
    plot_list(plt,glob.glob('smp_*data'),type='cycle')
    plt.legend(loc='lower center',ncols=3)
    plot_write('plot_smp_cycle',html)
# ----------------------------------------------------------------------

html_stop_row(html)

# ======================================================================

html_stop_table(html)
html_stop(html)
