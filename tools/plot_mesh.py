#!/usr/bin/python
import fileinput
import matplotlib
import matplotlib.pyplot as plt

# 10:101
#           a = 1
# 1 1        a = 2
#0 2        a = 4
# :
#1 0.5      t = 0.5
# 0          t = 0.25
#
def decode_block(block_name):
    prefix="B"
    dimension = block_name.count("_") + 1
    indices = block_name[len(prefix):].split("_")
    lower=[0,0,0]
    upper=[0,0,0]
    level=[0,0,0]
    index_axis = 0
    while index_axis < len(indices):
        index = indices[index_axis]
        index_array = 0
        index_tree  = 0
        c=index.find(":")
        index_forest = index[0:c] + index[c+1:]
        in_tree = 0
        if (index_forest==-1): in_tree = 1
        index_bit=0
        factor_array = 1
        factor_tree  = 1
        level = 0
        while index_bit < len(index):
            bit = index[index_bit]
            if (bit == ":"):
                in_tree = 1
                index_bit = index_bit + 1
            else:
                if (in_tree == 1): 
                    level = level + 1
                    factor_tree = 0.5*factor_tree
                    index_tree = index_tree + int(bit)*factor_tree
                else:
                    index_array = 2*index_array + int(bit)
                index_bit = index_bit + 1
          
        lower[index_axis] = index_array + index_tree
        upper[index_axis] = lower[index_axis] + factor_tree
        index_axis = index_axis + 1

    return [level,lower, upper]

fig = plt.figure()
ax = fig.add_subplot(111)
xmin=10000; xmax=0
ymin=10000; ymax=0
colormap = ['red','orange','green','blue','magenta','cyan']
lines = []
for line in fileinput.input():
      lines.append(line)
max_level = 10

for plotlevel in xrange(0,max_level):


      for line in lines:
        [level,lower, upper] = decode_block(line.split()[0])
        if (level == plotlevel):
            nc = len(colormap)
            r = matplotlib.patches.Rectangle ( (lower[0],lower[1]), 
                                               (upper[0]-lower[0]), 
                                               (upper[1]-lower[1]), 
                                               fill=False,
                                               edgecolor=colormap[level % nc])
            xmin = min(xmin,lower[0])
            ymin = min(ymin,lower[1])
            xmax = max(xmax,upper[0])
            ymax = max(ymax,upper[1])
            ax.add_patch(r)

xsize = xmax - xmin
ysize = ymax - ymin
dsize = abs(xsize - ysize)

# 1:1 aspect ratio
if (xsize > ysize):
    ymin -= 0.5*dsize
    ymax += 0.5*dsize
    ysize = ymax - ymin
elif (ysize > xsize):
    xmin -= 0.5*dsize
    xmax += 0.5*dsize
    xsize = xmax - xmin

# add border

r=0.2*xsize
xmin -= r
ymin -= r
xmax += r
ymax += r

plt.axes().set_aspect('equal','datalim')
    
plt.xlim([xmin,xmax])
plt.ylim([ymin,ymax])
plt.show()

