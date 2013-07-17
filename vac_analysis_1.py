# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from pymatgen.io.cifio import CifParser
from pymatgen.core.structure import Structure
from pymatgen.io.cifio import CifParser
from pymatgen.defects.point_defects import Vacancy
import pymatgen.defects.point_defects

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
for filen in filelist[0:20]:
    name = os.path.splitext(filen)[0]
    item = name.split('--')
    key = item[0]
    struct_descript_dict[key] = {'formula':''.join(item[1].split())}
    #print item

    #read the file to pymatgen struct
    filename = os.path.join(struct_dir,filen)
    try:
        struct = CifParser(filename).get_structures()[0]
    except:
        print filename
        del struct_descript_dict[key]
        continue
    try:
        vac = Vacancy(struct)
    except:
        print struct_descript_dict[key]['formula']
        del struct_descript_dict[key]
        continue
    no_vac = vac.defectsite_count()
    #print no_vac
    # add number of vacancies to descriptor
    struct_descript_dict[key]['no_vac'] = no_vac
    
    #Create a list storing a dictionary of vacancy properties for each vacancy
    vacancy_list = []
    valences = vac.valences.values()
    valences.sort()
    #print valences
    try:
        anion_cation_charge_ratio = abs(valences[-1]/valences[0])
    except:   #BVAnalyzer fails sometimes and valences are set to 0
        print struct_descript_dict[key]['formula']
        print valences
        anion_cation_charge_ratio = 0
    radii = vac.ionic_radii.values()
    radii.sort()
    try:
        anion_cation_radii_ratio = radii[-1]/radii[0]
    except:
        #print radii
        anion_cation_radii_ratio = 0
    for i in range(no_vac):
        vac_prop_dict = {}
        vac_prop_dict['coordination_no'] = vac.get_defectsite_coordination_number(i)
        vac_prop_dict['coordinated_elements'] = vac.get_coordinated_elements(i)
        vac_prop_dict['eff_charge'] = vac.get_defectsite_effective_charge(i)
        vac_prop_dict['el'] = vac.get_defectsite(i).species_string
        #vac_prop_dict['surf_area'] = vac.get_surface_area(i)
        #vac_prop_dict['volume'] = vac.get_volume(i)
        vac_prop_dict['chrg_ratio'] = anion_cation_charge_ratio
        vac_prop_dict['radii_ratio'] = anion_cation_radii_ratio
        vacancy_list.append(vac_prop_dict)
    struct_descript_dict[key]['vacancies'] = vacancy_list
        
#print struct_descript_dict


