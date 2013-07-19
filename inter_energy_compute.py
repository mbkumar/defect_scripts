#!/usr/bin/env python
# Generate the descriptors for intersitital defects

import os
import sys
import json

import matplotlib.pyplot as plt

from pymatgen.io.cifio import CifParser
from pymatgen.core.structure import Structure
from pymatgen.defects.point_defects import Interstitial, InterstitialFormationEnergy
from pymatgen.defects.point_defects import StructWithValenceIonicRadius
from pymatgen.symmetry.finder import SymmetryFinder


#cif_file = 'Binary oxide dataset/20130524/comp_structures/25779--Ti4 O6.cif'

def inter_energy_gen(cif_file):
    # Use the number as struct identifier
    # Store the remaining string as formula
    dir, filen = os.path.split(cif_file)
    name = os.path.splitext(filen)[0]
    item = name.split('--')
    #key = item[0]
    inter_descript_dict = {'formula':''.join(item[1].split())}

    #read the file to pymatgen struct
    try:
        struct = CifParser(cif_file).get_structures(False)[0]
    except:
        print 'reading ', cif_file, ' failed'
        return None, filen


    try:
        struct_val_rad = StructWithValenceIonicRadius(struct)
    except:
        print inter_descript_dict['formula']
        print "Unable to identify radii and valences"
        return None, filen

    inter = Interstitial(struct_val_rad)
    #inter.prune_defectsites('Mn', 3)
    no_inter = inter.defectsite_count()
    print "No of interstitials: ", no_inter
    #inter.reduce_defectsites()
    #no_inter = inter.defectsite_count()
    #print "No of interstitials: ", no_inter
    ife = InterstitialFormationEnergy(inter)
    inter_descript_dict['no_inter'] = no_inter

    #Create a list storing a dictionary of properties for each interstitial
    inter_list = []
    for i in range(no_inter):
        inter_dict = {}
        inter_dict['radius'] = inter.get_radius(i)
        #try:
        #    inter_dict['energy'] = ife.get_energy_fixsc(i, 'Mn', 3, relax=False)
        inter_dict['energy'] = ife.get_energy_fixsc(i, 'Mn', 3, relax=True)
        #except:
        #    inter_dict['energy'] = 'NaN'
        print i, inter_dict['radius'], inter_dict['energy']
        inter_list.append(inter_dict)
    inter_descript_dict['interstitials'] = inter_list
    inter_descript_dict

if __name__ == '__main__':
    Ti4O6_dict = inter_energy_gen(
            'Binary oxide dataset/20130524/comp_structures/33647--Mn16 O24.cif'
            )
    print Ti4O6_dict
    dump = json.dumps(Ti4O6_dict)
    with open('Ti4O6_energy.json', 'w') as fp:
        fp.write(dump)

