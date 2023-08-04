from nanocore import *
from nanocore.io import read_xsf
from nanocore.simulator import siesta as s2
import numpy as np
import copy
import os


os.chdir('input')
model = s2.read_struct("STRUCT.fdf")
sim = s2.Siesta(model)
sim.read_fdf("BASIS.fdf")
sim.read_fdf("RUN.fdf")
sim.read_fdf("KPT.fdf")
sim.read_fdf("TS.fdf")
os.chdir('..')


# read structures from external files
electrode = io.read_xsf('Au111.xsf') # 2d pbc on xz plane
channel = io.read_xsf('2H.xsf') # 3d pbc
supercell=(3,4)
d_x=(2, 2.5, 0.2)
d_y=(2, 2.5, 0.2)
d_z=(2, 2.5, 0.2)

# get cell matching
n1, n2 = electrode.get_cell_match(channel, 'y', 'y', max_unit=5)

# set repeat unit (USER INPUT)
m1, m2 = supercell[0], supercell[1]

# generate supercells
electrode_2 = electrode * (m1, n1, 1)
channel_2 = channel * (1, n2, m2)

# adjust cell size
ratio = channel_2.get_cell()[1][1] /electrode_2.get_cell()[1][1]
electrode_3 = electrode_2.adjust_cell_size(ratio, 7)

direction = ['x', 'y', 'z']
dir_range = [d_x, d_y, d_z]

for k in range(0,3):
    os.system('mkdir %s_dir'%direction[k])
    os.chdir('%s_dir'%direction[k])
    dir = dir_range[k]
    for d in np.arange(dir[0],dir[1],dir[2]):
        elec_left = copy.deepcopy(electrode_3)
        elec_right = copy.deepcopy(electrode_3)
        chan = copy.deepcopy(channel_2)
        a = 2
        b = 2
        c = 2
        if k == 0:
            a = d
        elif k == 1:
            b = d
        else:
            c = d

        os.system('mkdir %.2f'%d)
        zmax = elec_left.get_zmax()       # location of the surface
        chan.select_all()
        chan.translate(a, b, zmax+c) # move MoS2 on the surface
        zmax = chan.get_zmax()
        elec_right.select_all()
        elec_right.translate(0, 0, zmax+c)

        # add two structures
        new = elec_left + chan + elec_right  # merge two systems with PBC of "slab"

        # set vacuum along z axis
        zmax = elec_right.get_zmax()
        new.set_vacuum(zmax-new._cell[2][2]+15, 'z')

        # sort atoms by z coordinates
        new.select_all()
        new.sort('z')


        # save the interface model
        io.write_xsf('model.xsf', new)
        os.system('mv model.xsf %.2f'%d)

        # run simulator
        os.chdir('%.2f'%d)
        model = read_xsf("model.xsf")
        sim._atoms = model
        os.system('cp ~/pseudopotential/GGA/Mo.psf .')
        os.system('cp ~/pseudopotential/GGA/Au.psf .')
        os.system('cp ~/pseudopotential/GGA/S.psf .')
        sim.run("elec", 12)
        os.chdir('..')

    os.chdir('..')


