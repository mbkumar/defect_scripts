from math import sqrt
import json     #Get rid of json related function in future

_struct_prop = {'charge_ratio', 'radius_ratio', 'no_inter'}

def minkowski (struct1, struct2, power):
    """
    Computes the Minkowski distance (interstitial properties) between two structures
    """
    dist = 0
    for key in _struct_prop:
        dist += (abs(struct1[key] - struct2[key]))**power

    return dist**(1.0/power)


def manhattan (struct1, struct2):
    """
    Computes the Manhattan distance (interstitial properties) between two structures
    """
    dist = 0
    for key in _struct_prop:
        dist += abs(struct1[key] - struct2[key])

    return dist


def compute_nearest_neighbour(struct, structs):
    """
    Compute the Euclidean distance between the input structure with other 
    structures and sort the structures list with the closest one at the front.
    """
    distances = []
    for key in structs.keys():
        if key != struct:
            dist = minkowski(structs[struct], structs[key], 2.0)
            distances.append((distance, key))
            
    distances.sort()
    return distances


#Note: For structures, the different proprerties are completely
#unrelated. Instead it seems better to compute the pearsson coefficient
#between different properties
def pearson_coeff_struct(struct1, sturct2):
    """
    Compute the Pearsson coefficient between two structures.
    """
    dotprod = 0
    struct1_sum = 0
    struct2_sum = 0
    struct1_sq_sum = 0
    struct2_sq_sum = 0
    for key in _struct_prop:
        dotprod += struct1[key]*struct2[key]
        struct1_sum += struct1[key]
        struct2_sum += struct2[key]
        struct1_sq_sum += struct1[key]**2
        struct2_sq_sum += struct2[key]**2

    n = len(_struct_prop)
    r = dotprod - struct1_sum*struct2_sum/n
    r /= sqrt(struct1_sq_sum - (struct1_sum**2)/n)  
    r /= sqrt(struct1_sq_sum - (struct1_sum**2)/n)

    return r


def pearson_coeff1(prop1, prop2, structs_dict):
    """
    Compute the Pearsson coefficient between two properties of structures.
    """
    sd = structs_dict
    dotprod = 0
    prop1_sum = 0
    prop2_sum = 0
    prop1_sq_sum = 0
    prop2_sq_sum = 0
    for key in sd.keys():
        dotprod += sd[key][prop1]*sd[key][prop2]
        prop1_sum += sd[key][prop1]
        prop2_sum += sd[key][prop2]
        prop1_sq_sum += sd[key][prop1]**2
        prop2_sq_sum += sd[key][prop2]**2

    n = len(sd.keys())
    r = dotprod - prop1_sum*prop2_sum/n
    r /= sqrt(prop1_sq_sum - (prop1_sum**2)/n) 
    r /= sqrt(prop1_sq_sum - (prop1_sum**2)/n)

    return r

def pearson_coeff(prop1, prop2, structs_dict):
    """
    Compute the Pearsson coefficient between two properties of structures.
    """
    sd = structs_dict
    prop1_sum = 0
    prop2_sum = 0
    for key in sd.keys():
        prop1_sum += sd[key][prop1]
        prop2_sum += sd[key][prop2]
    n = len(sd.keys())
    prop1_mean = prop1_sum/n
    prop2_mean = prop2_sum/n
    dotprod = 0
    prop1_dlta_sq_sum = 0
    prop2_dlta_sq_sum = 0
    for key in sd.keys():
        prop1_dlta = sd[key][prop1] - prop1_mean
        prop2_dlta = sd[key][prop2] - prop2_mean
        dotprod += prop1_dlta*prop2_dlta
        prop1_dlta_sq_sum += prop1_dlta**2
        prop2_dlta_sq_sum += prop2_dlta**2

    r = dotprod / sqrt(prop1_dlta_sq_sum) / sqrt(prop2_dlta_sq_sum)

    return r

def pearson_coeff_arrays(X, Y):
    """
    Compute the Pearsson coefficient between two properties of structures.
    """
    X_sum = sum(X)
    Y_sum = sum(Y)
    assert(len(X) == len(Y))
    n = len(X)
    X_mean = X_sum/n
    Y_mean = Y_sum/n
    dotprod = 0
    X_dlta_sq_sum = 0
    Y_dlta_sq_sum = 0
    for i in range(n):
        x_dlta = X[i] - X_mean
        y_dlta = Y[i] - Y_mean
        dotprod += x_dlta*y_dlta
        X_dlta_sq_sum += x_dlta**2
        Y_dlta_sq_sum += y_dlta**2

    r = dotprod / sqrt(X_dlta_sq_sum) / sqrt(Y_dlta_sq_sum)

    return r


def pearson_coeff_onearray(X, Y):
    """
    Compute the Pearsson coefficient between two properties of structures.
    """
    X_sum = sum(X)
    Y_sum = sum(Y)
    assert(len(X) == len(Y))
    n = len(X)
    X_mean = X_sum/n
    Y_mean = Y_sum/n
    dotprod = 0
    X_dlta_sq_sum = 0
    Y_dlta_sq_sum = 0
    for i in range(n):
        x_dlta = X[i] - X_mean
        y_dlta = Y[i] - Y_mean
        dotprod += x_dlta*y_dlta
        X_dlta_sq_sum += x_dlta**2
        Y_dlta_sq_sum += y_dlta**2

    r = dotprod / sqrt(X_dlta_sq_sum) / sqrt(Y_dlta_sq_sum)

    return r


def return_dict(json_file):
    with open(json_file, 'r') as fp:
        json_string = fp.read()

    dictionary = json.loads(json_string)
    return dictionary


if __name__ == "__main__":
    inter_dict = return_dict('inter_prop.json')
    print json.dumps(inter_dict, sort_keys=False, indent=4)

