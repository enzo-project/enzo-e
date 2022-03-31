import yt
from yt.units import cm, g
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import use
use('Agg')

Nlist = [0,1]

@yt.derived_field(name="Etot", sampling_type="cell", units="erg")
def Etot(field, data):
    return data['total_energy'].to('erg/g')*data['density'].to('g/cm**3')*data['cell_volume'].to('cm**3') 
@yt.derived_field(name="metal_mass", sampling_type="gas", units="g")
def metal_tot(field, data):
    return data['density'].to('g/cm**3').d / data['density'].to('code_density').d * data['metal_density']*g*cm**-3*data['cell_volume'].to('cm**3')
Elist = []
Mlist = []
metalslist = []
for N in Nlist:
    ds = yt.load(f'FEEDBACK_TEST_{N:04d}/FEEDBACK_TEST_{N:04d}.block_list')
    ad = ds.all_data()
    Elist.append(ad.quantities.total_quantity('Etot'))
    Mlist.append(ad.quantities.total_quantity('cell_mass').to('Msun'))
    metalslist.append(ad.quantities.total_quantity('metal_mass').to('Msun'))

#print(Mlist)
print(f'Mass deposited: {Mlist[-1]-Mlist[0]}')
#print(metalslist)
print(f'Metals deposited: {metalslist[-1]-metalslist[0]}')
#print(Elist)
print(f'Energy deposited: {Elist[-1] - Elist[0]}')
