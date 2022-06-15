#import time
import numpy as np
import pandas as pd
from collections import defaultdict

import multiprocessing
from multiprocessing import Pool

# Custom library for performing spectroscopy
import pipeline_func

def tree():
    return defaultdict(tree)

#start_time = time.time()

# Location of directory containing FITS files
rootdir = "/Users/thepoetoftwilight/Documents/SDSS/CLIFs/"

from os import listdir
from os.path import isfile, join

# Load FITS files, reject any ones that are not
obj_list = [f for f in listdir(rootdir+"Spectra/") if isfile(join(rootdir+"Spectra/", f))]
obj_list = [f for f in obj_list if f[-5:]==".fits"]
obj_list = [f[0:-5] for f in obj_list]
n_obj = len(obj_list)


line_ids = ['H_alpha', 'H_beta', '[O_III] 5007', 'C_III] 1908', 'He_II 4685', '[Ne_III] 3868', '[N_II] 6583']
            
#line_ids = ['H_gamma', 'H_delta', 'H_epsilon', '[O_II] 3725', '[O_II] 3727', '[O_III] 4363', '[O_III] 4959']

# H-beta - flux, EW, cont.
emission_ew = tree()

#for obj_name in obj_list:
#    emission_ew[obj_name] = {}
#    for l in line_ids:
#        emission_ew[obj_name][l] = [0, 0]

#start_time = time.time()

if __name__ == "__main__":

    n_proc = multiprocessing.cpu_count()
    
    #print("About to start:")
    
    manager = multiprocessing.Manager()
    q = manager.Queue()

    with Pool(processes = n_proc) as pool:
        pool.map(pipeline_func.emission_pipeline, [(rootdir, obj, q, line_ids) for obj in obj_list])
        
    while not q.empty():
        output = q.get()
        obj_name = output[0]
        emission_dict = output[1]
        #print(obj_name + " done!")
        emission_ew[obj_name] = emission_dict

fluxes = tree()
        
ews = tree()

for l in line_ids:

    fluxes[l] = tree()
    ews[l] = tree()

    for obj_name in obj_list:

        fluxes[l][obj_name] = emission_ew[obj_name][l][0]
        ews[l][obj_name] = emission_ew[obj_name][l][1]

# Store emission line parameters in a .csv file
csv_dat_compile = []

for i in range(0, n_obj):

    obj_name = obj_list[i]

    csv_dat = [obj_name]

    for l in line_ids:
        csv_dat.append(fluxes[l][obj_name])
        csv_dat.append(ews[l][obj_name])

    csv_dat = np.array(csv_dat)
    
    csv_dat_compile.append(csv_dat)
    
csv_dat_compile = np.array(csv_dat_compile)

df = pd.DataFrame(csv_dat_compile)

cols = ['Object']

for l in line_ids:
    flux_label = l + ' Flux'
    ew_label = l + ' EW'
    cols.append(flux_label)
    cols.append(ew_label)

# Guarding for script runs
if(df.empty==False):
    df.columns = cols

df.to_csv(rootdir + "pipeline_emission.csv", index=False)

#df.to_csv(rootdir + "pipeline_emission_2.csv", index=False)

#print("--- %s seconds ---" % (time.time() - start_time))
