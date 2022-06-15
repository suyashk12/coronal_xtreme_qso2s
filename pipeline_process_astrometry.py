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
rootdir = "/Users/thepoetoftwilight/Documents/SDSS/Code/CLIFs/"

from os import listdir
from os.path import isfile, join

# Load FITS files, reject any ones that are not
obj_list = [f for f in listdir(rootdir+"Spectra/") if isfile(join(rootdir+"Spectra/", f))]
obj_list = [f for f in obj_list if f[-5:]==".fits"]
obj_list = [f[0:-5] for f in obj_list]
n_obj = len(obj_list)

# Object dictionary
obj_data = tree()

#start_time = time.time()

# Process all objects parallely
if __name__ == "__main__":

    n_proc = multiprocessing.cpu_count()

    #print("About to start:")

    manager = multiprocessing.Manager()
    q = manager.Queue()

    with Pool(processes = n_proc) as pool:
        pool.map(pipeline_func.process_obj_astrometry, [(rootdir, obj, q) for obj in obj_list])

    while not q.empty():
        output = q.get()
        obj_name = output[0]
        obj_dict = output[1]
        obj_data[obj_name] = obj_dict
        
    #print("All done!")

# Store emission line parameters in a .csv file
csv_dat_compile = []

for i in range(0, n_obj):

    obj_name = obj_list[i]

    csv_dat = []

    #Positional information
    csv_dat.append(obj_name)
    csv_dat.append(obj_data[obj_name]['ra'])
    csv_dat.append(obj_data[obj_name]['dec'])
    csv_dat.append(obj_data[obj_name]['redshift'])

    csv_dat = np.array(csv_dat)
    
    csv_dat_compile.append(csv_dat)
    
csv_dat_compile = np.array(csv_dat_compile)

df = pd.DataFrame(csv_dat_compile)

cols = ['Object', 'RA', 'Dec', 'z']

# Guarding for script runs
if(df.empty==False):
    df.columns = cols

df.to_csv(rootdir + "pipeline_astrometry.csv", index=False)

#df.to_csv(rootdir + "pipeline_emission_2.csv", index=False)

#print("--- %s seconds ---" % (time.time() - start_time))
