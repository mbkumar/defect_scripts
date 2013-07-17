#!/usr/bin/env python
# Generate the descriptors for intersitital defects

import os
import sys
import json

import matplotlib.pyplot as plt
from pymongo import MongoClient

from pymatgen.io.cifio import CifParser
from pymatgen.core.structure import Structure
from pymatgen.defects.point_defects import Interstitial
from pymatgen.defects.point_defects import StructWithValenceIonicRadius
from pymatgen.symmetry.finder import SymmetryFinder


#cif_dir = 'Binary oxide dataset/20130524/comp_structures'
def gen_inter_descript(cif_dir, collection, prune_args):
    """
    Generates the interstitial descriptor values for the structures 
    in input cif directory and store them into the database supplied
    Args:
        cif_dir:
            Directory containing the cif files 
        collection:
            Mongo database collection
        prune_ars:
            List of keywords to prune the interstitials
            a) anion:
                Only the interstitial sites that can be occupied by the 
                anion are kept
            a) cation:
                Only the interstitial sites that can be occupied by the 
                smallest cation are kept
            c) [El, oxi_state]:
                Only the interstitial sites that can be occupied by the 
                input ion are kept
                Ex: ['Li', 1]
    """
    filelist = os.listdir(cif_dir)
    for filen in filelist:
        if os.path.splitext(filen)[1] != '.cif':
            filelist.remove(filen)

    for filen in filelist:
        if os.path.splitext(filen)[1] != '.cif':
            filelist.remove(filen)

    # Use the number as struct identifier
    # Store the remaining string as formula
    failed_struct = []
    for filen in filelist:
        filename = os.path.join(cif_dir,filen)
        status = inter_descript_gen(filename, collection, prune_args)
        if not status:
            failed_struct.append(filen)

    return failed_struct

def inter_descript_gen(cif_file, collection, prune_args):
    """
    Generates the interstitial descriptor values for the structures 
    in input cif directory and store them into the database supplied
    Args:
        cif_file:
            cif File 
        collection:
            Mongo database collection
        prune_ars:
            List of keywords to prune the interstitials
            a) anion:
                Only the interstitial sites that can be occupied by the 
                anion are kept
            a) cation:
                Only the interstitial sites that can be occupied by the 
                smallest cation are kept
            c) [El, oxi_state]:
                Only the interstitial sites that can be occupied by the 
                input ion are kept
                Ex: ['Li', 1]
    """
    # Use the number as struct identifier
    # Store the remaining string as formula
    print cif_file
    dir, filen = os.path.split(cif_file)
    print filen
    name = os.path.splitext(filen)[0]
    item = name.split('--')
    key = item[0]
    inter_descript_dict = {'key': key, 'formula':''.join(item[1].split())}
    FAIL=0
    SUCC=1

    #read the file to pymatgen struct
    try:
        struct = CifParser(cif_file).get_structures(False)[0]
    except:
        print 'reading ', cif_file, ' failed'
        return FAIL

    symm_finder = SymmetryFinder(struct)
    inter_descript_dict['structure'] = struct
    inter_descript_dict['crystal_system'] = symm_finder.get_crystal_system()
    inter_descript_dict['spacegroup_no'] = symm_finder.get_spacegroup_number()
    no_symmops = len(symm_finder.get_symmetry_operations())
    inter_descript_dict['no_symmops'] = no_symmops

    try:
        struct_val_rad = StructWithValenceIonicRadius(struct)
        valences = struct_val_rad.valences.values()
        valences.sort()
        anion_cation_charge_ratio = abs(valences[-1]/valences[0])
        radii = struct_val_rad.radii.values()
        radii.sort()
        anion_cation_radii_ratio = radii[-1]/radii[0]
    except:
        print inter_descript_dict['formula']
        print "Unable to identify radii and valences"
        return FAIL
    inter_descript_dict['charge_ratio'] = anion_cation_charge_ratio 
    inter_descript_dict['radius_ratio'] = anion_cation_radii_ratio 

    try:
        inter = Interstitial(struct_val_rad)
    except:
        print inter_descript_dict['formula']
        print 'interstitial not generated'
        return FAIL


    if prune_args[0] == "cation":
        inter.radius_prune_defectsites(radii[0])    #smallest element in the structure
    elif prune_args[0] == "anion":
        inter.radius_prune_defectsites(radii[-1])   #largest element in the structure
    else:
        el = prune_args[0]
        oxi_state = prune_args[1]
        inter.prune_defectsites(el, oxi_state)
    inter.prune_close_defectsites()

    if inter.defectsite_count() > 10:
        inter.reduce_defectsites()
    no_inter = inter.defectsite_count()
    inter_descript_dict['no_inter'] = no_inter
    #print >>sys.stderr, inter_descript_dict[key]['formula'], "No. of inter", no_inter

    #Create a list storing a dictionary of properties for each interstitial
    inter_list = []
    for i in range(no_inter):
        inter_dict = {}
        inter_dict['coord_no'] = inter.get_defectsite_coordination_number(i)
        inter_dict['coord_el'] = list(inter.get_coordinated_elements(i))
        inter_dict['radius'] = inter.get_radius(i)
        inter_dict['coord_chrg_sum'] = inter.get_coordsites_charge_sum(i)
        inter_list.append(inter_dict)
    inter_descript_dict['interstitials'] = inter_list
    descriptor_id = collection.insert(inter_descript_dict)
    return descriptor_id

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
    client = MongoClient()
    db = client['BinOxide_testbench_interstitials']

    #cation_inter_collection = db['cation_interstitial']
    #failed_struct = gen_inter_descript(
    #        'Binary oxide dataset/20130524/comp_structures', 
    #        cation_inter_collection,
    #        'cation'
    #        )

    #anion_inter_collection = db['anion_interstitial']
    #failed_struct = gen_inter_descript(
    #        'Binary oxide dataset/20130524/comp_structures', 
    #        anion_inter_collection,
    #        'anion'
    #        )
    #dump = json.dumps(inter_anion_prune_dict)
    #with open('inter_prop_anion_prune.json', 'w') as fp:
    #    fp.write(dump)

    Li_inter_collection = db['Li_interstitial']
    failed_struct = gen_inter_descript(
            'Binary oxide dataset/20130524/comp_structures', 
            Li_inter_collection,
            ('Li',1)
            )
