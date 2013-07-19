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


#cif_dir = 'Binary oxide dataset/20130524/comp_structures'
def inter_descript_generator(cif_dir):
    filelist = os.listdir(cif_dir)
    for filen in filelist:
        if os.path.splitext(filen)[1] != '.cif':
            filelist.remove(filen)

    for filen in filelist:
        if os.path.splitext(filen)[1] != '.cif':
            filelist.remove(filen)

    # Use the number as struct identifier
    # Store the remaining string as formula
    inter_descript_dict = {}
    failed_struct = []
    for filen in filelist:
        name = os.path.splitext(filen)[0]
        item = name.split('--')
        key = item[0]
        inter_descript_dict[key] = {'formula':''.join(item[1].split())}

        #read the file to pymatgen struct
        filename = os.path.join(cif_dir,filen)
        try:
            struct = CifParser(filename).get_structures(False)[0]
        except:
            print 'reading ', filename, ' failed'
            del inter_descript_dict[key]
            failed_struct.append(filen)
            continue

        try:
            struct_val_rad = StructWithValenceIonicRadius(struct)
        except:
            print inter_descript_dict[key]['formula']
            print "Unable to identify radii and valences"
            del inter_descript_dict[key]
            failed_struct.append(filen)
            continue

        try:
            inter = Interstitial(struct_val_rad)
        except:
            print inter_descript_dict[key]['formula']
            print 'interstitial not generated'
            del inter_descript_dict[key]
            failed_struct.append(filen)
            continue
        inter.prune_defectsites()
        no_inter = inter.defectsite_count()
        inter_descript_dict[key]['no_inter'] = no_inter
        #print >>sys.stderr, inter_descript_dict[key]['formula'], "No. of inter", no_inter
        
        #Create a list storing a dictionary of properties for each interstitial
        valences = inter.struct_valences.values()
        valences.sort()
        try:
            anion_cation_charge_ratio = abs(valences[-1]/valences[0])
        except:   #BVAnalyzer fails sometimes and valences are set to 0
            print inter_descript_dict[key]['formula']
            print valences
            anion_cation_charge_ratio = 0
        inter_descript_dict[key]['charge_ratio'] = anion_cation_charge_ratio 

        radii = inter.struct_radii.values()
        radii.sort()
        try:
            anion_cation_radii_ratio = radii[-1]/radii[0]
        except:
            anion_cation_radii_ratio = 0
        inter_descript_dict[key]['radius_ratio'] = anion_cation_radii_ratio 

        inter_list = []
        for i in range(no_inter):
            inter_dict = {}
            inter_dict['coord_no'] = inter.get_defectsite_coordination_number(i)
            inter_dict['coord_el'] = list(inter.get_coordinated_elements(i))
            inter_dict['radius'] = inter.get_radius(i)
            inter_dict['coord_chrg_sum'] = inter.get_coordsites_charge_sum(i)
            inter_list.append(inter_dict)
        inter_descript_dict[key]['interstitials'] = inter_list
    return inter_descript_dict, failed_struct
        
def gen_inter_descript(cif_dir, prune_args):
    filelist = os.listdir(cif_dir)
    for filen in filelist:
        if os.path.splitext(filen)[1] != '.cif':
            filelist.remove(filen)

    for filen in filelist:
        if os.path.splitext(filen)[1] != '.cif':
            filelist.remove(filen)

    # Use the number as struct identifier
    # Store the remaining string as formula
    inter_descript_dict = {}
    failed_struct = []
    for filen in filelist:
        filename = os.path.join(cif_dir,filen)
        key, inter_dict = inter_descript_gen(filename, prune_args)
        if key:
            inter_descript_dict[key] = inter_dict
        else:
            failed_struct.append(filen)

    return inter_descript_dict, failed_struct

def inter_descript_gen(cif_file, prune_args):
    # Use the number as struct identifier
    # Store the remaining string as formula
    print cif_file
    dir, filen = os.path.split(cif_file)
    print filen
    name = os.path.splitext(filen)[0]
    item = name.split('--')
    key = item[0]
    inter_descript_dict = {'formula':''.join(item[1].split())}

    #read the file to pymatgen struct
    try:
        struct = CifParser(cif_file).get_structures(False)[0]
    except:
        print 'reading ', cif_file, ' failed'
        return None, filen

    symm_finder = SymmetryFinder(struct, symprec=1e-1)
    inter_descript_dict['crystal_system'] = symm_finder.get_crystal_system()
    inter_descript_dict['spacegroup_no'] = symm_finder.get_spacegroup_number()
    no_symmops = len(symm_finder.get_symmetry_operations())
    inter_descript_dict['no_symmops'] = no_symmops

    try:
        struct_val_rad = StructWithValenceIonicRadius(struct)
    except:
        print inter_descript_dict['formula']
        print "Unable to identify radii and valences"
        return None, filen

    try:
        inter = Interstitial(struct_val_rad)
    except:
        print inter_descript_dict['formula']
        print 'interstitial not generated'
        return None, filen


    
    valences = inter.struct_valences.values()
    valences.sort()
    try:
        anion_cation_charge_ratio = abs(valences[-1]/valences[0])
    except:   #BVAnalyzer fails sometimes and valences are set to 0
        print inter_descript_dict['formula']
        print valences
        anion_cation_charge_ratio = 0
    inter_descript_dict['charge_ratio'] = anion_cation_charge_ratio 

    radii = inter.struct_radii.values()
    radii.sort()
    try:
        anion_cation_radii_ratio = radii[-1]/radii[0]
    except:
        anion_cation_radii_ratio = 0
    inter_descript_dict['radius_ratio'] = anion_cation_radii_ratio 


    if prune_args[0] == "cation":
        inter.radius_prune_defectsites(radii[0])    #smallest element in the structure
    elif prune_args[0] == "anion":
        inter.radius_prune_defectsites(radii[-1])   #largest element in the structure
    else:
        el = prune_args[0]
        oxi_state = prune_args[1]
        print el, oxi_state
        inter.prune_defectsites(el, oxi_state)
    inter.prune_close_defectsites()

    if inter.defectsite_count() > 10:
    # Something wrong with the symmetry reduction of defects. 
    # Ignoring the outlier
        return None, filen
        inter.reduce_defectsites()
    no_inter = inter.defectsite_count()
    if no_inter > 5:
        print "Tag:", key, inter_descript_dict['formula'], no_inter
    inter_descript_dict['no_inter'] = no_inter
    #print >>sys.stderr, inter_descript_dict[key]['formula'], "No. of inter", no_inter
    #ife = InterstitialFormationEnergy(inter)

    #Create a list storing a dictionary of properties for each interstitial
    inter_list = []
    for i in range(no_inter):
        inter_dict = {}
        inter_dict['coords'] = inter.get_defectsite(i).coords
        inter_dict['coord_no'] = inter.get_defectsite_coordination_number(i)
        inter_dict['coord_el'] = list(inter.get_coordinated_elements(i))
        inter_dict['radius'] = inter.get_radius(i)
        inter_dict['coord_chrg_sum'] = inter.get_coordsites_charge_sum(i)
        inter_list.append(inter_dict)
    inter_descript_dict['interstitials'] = inter_list
    return key, inter_descript_dict

# define different plot generators
def coord_no_vs_inter_radius(inter_descript_dict):
    X = []
    Y = []
    dct = inter_descript_dict
    for key in dct.keys():
        for i in range(dct[key]['no_inter']):
            x = dct[key]['interstitials'][i]['coord_no'] 
            y = dct[key]['interstitials'][i]['radius']
            X.append(x)
            Y.append(y)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(X,Y)
    plt.savefig('test.png')

if __name__ == '__main__':
    inter_cation_prune_dict, failed_struct = gen_inter_descript(
            'Binary oxide dataset/20130524/comp_structures', ('cation',)
            )
    dump = json.dumps(inter_cation_prune_dict)
    with open('inter_prop_cation_prune.json', 'w') as fp:
        fp.write(dump)

    #inter_anion_prune_dict, failed_struct = gen_inter_descript(
    #        'Binary oxide dataset/20130524/comp_structures', ('anion',)
    #        )
    #dump = json.dumps(inter_anion_prune_dict)
    #with open('inter_prop_anion_prune.json', 'w') as fp:
    #    fp.write(dump)

    #inter_Li_prune_dict, failed_struct = gen_inter_descript(
    #        'Binary oxide dataset/20130524/comp_structures', ('Li',1)
    #        )
    #dump = json.dumps(inter_Li_prune_dict)
    #with open('inter_prop_Li_prune.json', 'w') as fp:
    #    fp.write(dump)
    print failed_struct
    #print struct_inter_dict
    #coord_no_vs_inter_radius(struct_inter_dict)
