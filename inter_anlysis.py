
import os
import sys
import json

import matplotlib.pyplot as plt

from pymatgen.io.cifio import CifParser
from pymatgen.core.structure import Structure
from pymatgen.defects.point_defects import Interstitial

#cif_dir = 'Binary oxide dataset/20130524/comp_structures'
        
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
    plt.savefig('test1.png')

if __name__ == '__main__':
    with open('inter_prop.json', 'r') as fp:
        json_string = fp.read()
    inter_descript_dict = json.loads(json_string)
    coord_no_vs_inter_radius(struct_inter_dict)
