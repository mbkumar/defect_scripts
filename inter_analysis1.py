# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from pymatgen.io.cifio import CifParser
from pymatgen.core.structure import Structure
from pymatgen.defects.point_defects import Interstitial

struct_dir = 'Binary oxide dataset/20130524/comp_structures'
filelist = os.listdir(struct_dir)
for filen in filelist:
    if os.path.splitext(filen)[1] != '.cif':
        filelist.remove(filen)

for filen in filelist:
    if os.path.splitext(filen)[1] != '.cif':
        filelist.remove(filen)

# Use the number as struct identifier
# Store the remaining string as formula
struct_descript_dict = {}
for filen in filelist[0:10]:
    name = os.path.splitext(filen)[0]
    item = name.split('--')
    key = item[0]
    struct_descript_dict[key] = {'formula':''.join(item[1].split())}
    print name

    #read the file to pymatgen struct
    filename = os.path.join(struct_dir,filen)
    try:
        struct = CifParser(filename).get_structures(primitive=False)[0]
    except:
        print 'reading ', filename, ' failed'
        del struct_descript_dict[key]
        print struct_descript_dict.keys()
        continue

    import sys
    #print >>sys.stderr, struct
    try:
        inter = Interstitial(struct)
    except:
        print 'interstitial not generated for ', struct_descript_dict[key]['formula']
        del struct_descript_dict[key]
        continue

