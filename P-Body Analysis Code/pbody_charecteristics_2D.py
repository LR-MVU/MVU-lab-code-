#################
#Import Packages
#################
import numpy as np
import pandas as pd
import os
from bigfishUI import *
from tifffile import imread
from skimage import measure
from skimage.measure import regionprops
import csv
from math import pi


##########
#Functions
##########
def get_filenames(directory, chan, filetype = '.tif'):
    '''gets all filenames corresponding to input channel
    str, str -> 2D list of match objects'''
    
    tif_pattern = re.compile(".*" + chan + ".*" + filetype)
    tif_files = [f for f in os.listdir(directory) if re.match(tif_pattern, f)]

    tif_files.sort()
    return tif_files

def write_csv(file_name, array, headers):
    with open(file_name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(array)
        
def format_results(arr):
    # If all conditions have no data
    if not any(len(row) for row in arr):
        return np.empty((0, len(arr)), dtype=object)
    
    max_len = max((len(row)) for row in arr)

    arr = np.array([
        list(row) + [None] * (max_len - (len(row)))
        for row in arr
    ], dtype=object) 
    return arr.T
    

def PBody_Charecteristics(conditions, diffraction, path, pixel_res, pbody_chan, cell_mask_chan = "print", nuc_mask_chan = "nuc", PB_folder = "Pb_Masks"):

    PBsize_dict = {} #dictionary to hold all PB sizes
    PBnum = {}
    PB_R_dict = {} #dictionary to hold all PB radii
    
    #Loop through all conditions
    for cond in conditions:
        print(cond)
        path_input = os.path.join(path, cond)
        path_PB = os.path.join(path_input, PB_folder)
        #get lists of files
        cell_mask_file = get_filenames(path_input, cell_mask_chan)
        nuc_mask_file = get_filenames(path_input, nuc_mask_chan)
        PB_mask_file = get_filenames(path_PB, pbody_chan)
        
        #Inialize lists to store data for the condition
        PB_areas = [] # Area in um^2
        PB_Rs = [] # Radius in um
        PB_num = [] # Number of PBodies

        #Loop through each cell
        for i in range(len(PB_mask_file)):
            print(PB_mask_file[i])
            #Load data
            PB_mask = np.squeeze(imread(os.path.join(path_PB, PB_mask_file[i]))) > 0
            cell_mask = np.squeeze(imread(os.path.join(path_input, cell_mask_file[i]))) > 0
            nuc_mask = np.squeeze(imread(os.path.join(path_input, nuc_mask_file[i]))) > 0

            cyt_mask = cell_mask & (~nuc_mask)
            PB_mask = PB_mask & cyt_mask
            labeled_pb_mask = measure.label(PB_mask, connectivity=1)
            pb_props = regionprops(labeled_pb_mask)
            
            pb_num = 0
            #Loop through each PB
            for prop in pb_props:
                #Compute radius of PB
                radius = (prop.equivalent_diameter_area)*pixel_res/2 - diffraction
                if radius > 0:
                    PB_Rs.append(radius)
                    PB_areas.append(radius**2*pi)
                    pb_num += 1
            
            PB_num.append(pb_num)
    
        #Add data to dictionary with condition as key
        PBsize_dict[cond] = PB_areas
        PB_R_dict[cond] = PB_Rs
        PBnum[cond] = PB_num
    
    PB_size_results = []
    PB_R_results = []
    PBnum_results = []
    for cond in conditions:
        PB_size_results.append(PBsize_dict[cond])
        PB_R_results.append(PB_R_dict[cond])
        PBnum_results.append(PBnum[cond])
    
    f_PB_size_results = format_results(PB_size_results)
    f_PB_R_results = format_results(PB_R_results)
    f_PBnum_results = format_results(PBnum_results)
    
    df_PB_size_results = pd.DataFrame(f_PB_size_results, columns=conditions)
    df_PB_R_results = pd.DataFrame(f_PB_R_results, columns=conditions)
    df_PBnum_results = pd.DataFrame(f_PBnum_results, columns=conditions)
    
    output_path = os.path.join(path, "Pbody Charecteristics " + pbody_chan +".xlsx")
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_PB_size_results.to_excel(writer, sheet_name="PB Size", index=False)
        df_PB_R_results.to_excel(writer, sheet_name="PB Radius", index=False)
        df_PBnum_results.to_excel(writer, sheet_name="PB Number", index=False)
    
    













