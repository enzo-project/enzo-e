import os
import sys
import h5py

'''
Need to put root blocks first
'''
name = 'ENZOE_LYAF-0120'
prefix = 'writer'
max_level = 4

outfile = open(f'{name}/{name}.block_list','w')

for level in range(max_level+1):
    for fname in os.listdir(name):
        #counter += 1
        if (prefix not in fname):
            continue

        #print(f'{name}/{fname}; {counter}/{len(os.listdir(name))}')

        for l in h5py.File(f'{name}/{fname}','r').keys():
            if level == 0 and (':' not in l): 
                outfile.write(f'{l} {fname}\n')
            elif level > 0 and (':' in l):
                if len(l[ l.rindex(':')+1: ]) == level:
                    outfile.write(f'{l} {fname}\n')     
                
outfile.close()
