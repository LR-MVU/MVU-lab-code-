# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:51:00 2024
@author: Lydia Hodgins

This code uses the Big FISH package to detect the position of both mRNA and P-Bodies from smFISH-IF images.

Before running make sure to change the channel variables to match the your inputed data (lines 51-53).
Also make sure the input and output folder names are correct (lines 46-47).
The output folder will be automatically created if it does not already exist.

For this code to run properly the images for each condition need to be in separate folders.
The names of the folder for each condition is defined in line 47 and assigned to the Condition variable.

Order of operations:
    1) All control P-Body images are loaded. The spot detection threshold is determined using all of these images. 
    The P-Body spots are detected for ctl.
    2) For all other conditions the P-Bodies are detected using the threshold determined from the ctl images.
    3) For each individual image the mRNA spot detection threshold is uniquely computed and the spots are detected.
    4) Dense regions of the mRNA spot detections are declustered.
    5) Cell mask and nucleus mask are loaded. Only spots that are within the cell and outside the nucleus are saved.
    6) P-Body and mRNA positions are saved in separate csv files. (One file for each image/cell)

"""

import numpy as np
import os
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.multistack as multistack
import bigfish.segmentation as segmentation
import bigfish.plot as plot
from matplotlib import pyplot as plt
import pandas as pd
import glob
import re
from bigfishUI import *
import csv

###### Functions #######

def write_csv(file_name, array, headers):
    with open(file_name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(array)
        
def csv_reader (csv_path):
   df_csv = pd.read_csv((csv_path), dtype = str)
   return df_csv.to_numpy()

def get_filenames(directory, chan):
    '''gets all filenames corresponding to input channel
    str, str -> 2D list of match objects'''
    
    lis = os.listdir(directory)
    
    tif_pattern = re.compile(".*_" + chan + ".*.tif")
    matches = regexMatchFilter(lis, tif_pattern)
    tif_files = [match.group() for match in matches]
    if chan == "555S":
        tif_files = [file for file in tif_files if not ("nuc" in file or "print" in file)]
    #print(tif_files)
    tif_files.sort()
    #print(tif_files)
    return tif_files 


##### CHANGE THESE VARIABLES #####
manual = False # SET TO FALSE FIRST RUN, CHANGE THRESHOLDS IF NEEDED IN SPREADSHEET, SET TO TRUE AND RUN AGAIN
run_all = True # Must be true for automatic and can be true if manual but running all again
manual_threshold_path = r"Z:\Morgane\MV198-199-200-201-202_25-01-22_MEFs HSP90AA Cy3 HSP110 Cy5 Dcp1A Cy7\Export\bigFISH_In\MEFs Hsp90aa - Hsp110_Thresholds.csv" # SET TO PATH TO SPREASHEET CREATED IN YOUR FOLDER
reference_threshold_path = r"" # Need a reference path if running manual and not running all
##### ********** #######

manual_counter = -1
if manual:
    raw = csv_reader(manual_threshold_path)
    input_thresholds = (raw[:,1:]).astype("float64")
if not run_all:
    raw = csv_reader(reference_threshold_path)
    reference_threshold = (raw[:,1:]).astype("float64")
else:
    reference_threshold = []
    
manual_thresholds = []
path_root = getImagesDir() #get path to directory containing cropped tif stacks
Conditions = [
               # "DMSO_CT",
               # "DMSO_HS",
               # "DMSO_2HR",
               # "DMSO_3HR",
               # "DMSO_4HR",
               # "DMSO_6HR",
               # "DMSO_8HR",
               # "HRN_CT",
               # "HRN_HS",
               # "HRN_2HR",
               # "HRN_3HR",
               # "HRN_4HR",
               # "HRN_6HR",
               # "HRN_8HR",
               # "Puro_CT",
               # "Puro_HS",
               # "Puro_2HR",
               # "Puro_3HR",
               # "Puro_4HR",
               # "Puro_6HR",
               # "Puro_8HR", 
               # "CT+ActD_CT",
               # "CT+ActD_3HR",
               # "CT+ActD_6HR",
               # "CT+ActD_9HR", 
               # "CT+ActD_12HR",
               # "ActD+HS_HS",
               # "ActD+HS_3HR",
               # "ActD+HS_6HR", 
               # "ActD+HS_9HR", 
               # "ActD+HS_12HR",
               # "HS+ActD_3HR",
               # "HS+ActD_6HR",
               # "HS+ActD_9HR",
               # "HS+ActD_12HR",
               # "HS+ActD_HS",
               #"Ctl",
               #"HS",
               #"HS_2H_R",
               #"HS_3H_R",
               #"HS_4H_R",
               #"HS_6H_R",
               #"HS_8H_R",
               #"HS_10H_R",
               #"HS_12H_R",
               #"HS_24H_R",
               "2hrs"
               ]
               
for cond in Conditions:
    name_list = []
    print(cond)
    path_input = os.path.join(path_root, cond)
    path_output = os.path.join(path_input, "output")
    
    if not os.path.exists(path_output): #make output folder if not already exists
        os.makedirs(path_output)
    
    DAPI = '395S'
    RNA_1 = '555S'
    #RNA_2 = '640S'
    IF = '640S' # CHANGE TO CORRECT CHANNEL
    cell_mask = 'print'
    nuc_mask = 'nuc'
    
    DAPI_files = get_filenames(path_input, DAPI)
    IF_files = get_filenames(path_input, IF)
    RNA_1_files = get_filenames(path_input, RNA_1)
    #RNA_2_files = get_filenames(path_input, RNA_2)
    cell_mask_files = get_filenames(path_input, cell_mask)
    nuc_mask_files = get_filenames(path_input, nuc_mask)
    
    threshold_file = os.path.join(path_output, "thresholds.csv")
    
    ##### Parameter Definitions #######
    
    # DECLUSTERING PARAMETERS
    alpha=0.5  # alpha impacts the number of spots per candidate region
    beta=1 # beta impacts the number of candidate regions to decompose
    gamma=7  # gamma the filtering step to denoise the image
    
    PB_spot_radius = (350, 270, 270) # in nanometer (one value per dimension zyx)
    
    RNA_spot_radius = (350, 150, 150) # in nanometer (one value per dimension zyx)
    
    Voxel_size = (300, 107.5, 107.5) # in nanometer (one value per dimension zyx)
    
    IF_imgs = []
    IF_mips = []

   
    for i in range(len(IF_files)):
        
        TIF_path = os.path.join(path_input, IF_files[i])
        IFfilename = IF_files[i][:-4]
        name_list.append(IFfilename)
        
        IF_im = stack.read_image(TIF_path)
        IF_imgs.append(IF_im)
        
        IF_mip = stack.maximum_projection(IF_im)
        IF_mips.append(IF_mip)
        
        
        
    #Detect PBodies using threshold for all images in control condition

    if not manual:

        if cond == "2hrs":
        
            plot.plot_elbow(
                images=IF_imgs, 
                voxel_size=Voxel_size, 
                spot_radius=PB_spot_radius,
                path_output=os.path.join(path_output, IFfilename+"_all_elbow"))  
            IF_spots, IF_threshold = detection.detect_spots(images=IF_imgs, return_threshold=True, 
                                                            voxel_size=Voxel_size,  
                                                            spot_radius=PB_spot_radius)   
            print(IF_threshold)
        
        else:
            IF_spots, IF_threshold = detection.detect_spots(images=IF_imgs, threshold = IF_threshold, return_threshold=True,
                                                                voxel_size=Voxel_size,  
                                                                spot_radius=PB_spot_radius)
            print(IF_threshold)
    
    
    for i in range(len(IF_files)):
        current_thresholds = []
        current_thresholds.append(name_list[i])
        
        manual_counter += 1
        
        if (not run_all) and manual and (input_thresholds[manual_counter] == reference_threshold[manual_counter]).all():
            continue
        
        TIF_1_path = os.path.join(path_input, RNA_1_files[i])
        RNA_1_filename = RNA_1_files[i][:-4]

        #TIF_2_path = os.path.join(path_input, RNA_2_files[i])
        #RNA_2_filename = RNA_2_files[i][:-4]
        
        RNA_1_im = stack.read_image(TIF_1_path)
        RNA_1_mip = stack.maximum_projection(RNA_1_im)
      
        #RNA_2_im = stack.read_image(TIF_2_path)
        #RNA_2_mip = stack.maximum_projection(RNA_2_im)        
        
        IF_im = IF_imgs[i]
        IF_mip = IF_mips[i]
        IFfilename = IF_files[i][:-4]
        
               
        #Detect mRNA using threshold for individual image
        plot.plot_elbow(
            images=RNA_1_im, 
            voxel_size=Voxel_size, 
            spot_radius=RNA_spot_radius,
            path_output=os.path.join(path_output, RNA_1_filename+"_elbow"))  
 
        # plot.plot_elbow(
        #     images=RNA_2_im, 
        #     voxel_size=Voxel_size, 
        #     spot_radius=RNA_spot_radius,
        #     path_output=os.path.join(path_output, RNA_2_filename+"_elbow"))

        if manual:
            if run_all or input_thresholds[manual_counter, 0] != reference_threshold[manual_counter, 0]:
                threshold_IF = input_thresholds[manual_counter, 0]
                
                plot.plot_elbow(
                    images=IF_im, 
                    voxel_size=Voxel_size, 
                    spot_radius=PB_spot_radius,
                    path_output=os.path.join(path_output, IFfilename+"_all_elbow"))
                
                IF_spot, IF_threshold = detection.detect_spots(images=IF_im, threshold = threshold_IF, return_threshold=True,
                                                                    voxel_size=Voxel_size,  
                                                                    spot_radius=PB_spot_radius)
            else:
                IF_threshold = input_thresholds[manual_counter, 0]
        else:
            IF_spot = IF_spots[i]
              
            
        current_thresholds.append(IF_threshold)    
            
        if not manual:
            RNA_1_spot, RNA_1_threshold = detection.detect_spots(images=RNA_1_im, return_threshold=True, 
                                                                 voxel_size=Voxel_size,  
                                                                 spot_radius=RNA_spot_radius)
        else:
            if run_all or input_thresholds[manual_counter, 1] != reference_threshold[manual_counter, 1]:
                threshold_RNA1 = input_thresholds[manual_counter, 1]
                RNA_1_spot, RNA_1_threshold = detection.detect_spots(images=RNA_1_im, return_threshold=True, threshold = threshold_RNA1,
                                                                    voxel_size=Voxel_size,  
                                                                    spot_radius=RNA_spot_radius)
            else: 
                RNA_1_threshold = input_thresholds[manual_counter, 1]
        current_thresholds.append(RNA_1_threshold)
        
        # if not manual:
        #     RNA_2_spot, RNA_2_threshold = detection.detect_spots(images=RNA_2_im, return_threshold=True, 
        #                                                          voxel_size=Voxel_size,  
        #                                                          spot_radius=RNA_spot_radius)
        # else:
        #     if run_all or input_thresholds[manual_counter, 2] != reference_threshold[manual_counter, 2]:
        #         threshold_RNA2 = input_thresholds[manual_counter, 2]
        #         RNA_2_spot, RNA_2_threshold = detection.detect_spots(images=RNA_2_im, return_threshold=True, threshold = threshold_RNA2,
        #                                                             voxel_size=Voxel_size,  
        #                                                             spot_radius=RNA_spot_radius)
        #     else: 
        #         RNA_2_threshold = input_thresholds[manual_counter, 2]
        # current_thresholds.append(RNA_2_threshold)
        # manual_thresholds.append(current_thresholds)
          
        
        ###### Get Nucleus Mask #######
        TIF_path = os.path.join(path_input, nuc_mask_files[i])
        filename = nuc_mask_files[i][:-4]
        mask_im = stack.read_image(TIF_path)
        print(TIF_path)
        print(mask_im.shape)
        mask_im = mask_im.astype(bool)
        nuc_label = segmentation.label_instances(mask_im)
        
        ###### Get Cell Masks ######
        TIF_path = os.path.join(path_input, cell_mask_files[i])
        filename = cell_mask_files[i][:-4]
        mask_im = stack.read_image(TIF_path)
        print(TIF_path)
        print(mask_im.shape)
        mask_im = mask_im.astype(bool)
        cell_label = segmentation.label_instances(mask_im)
        
        # Match nucleus mask with cell mask
        nuc_label, cell_label = multistack.match_nuc_cell(nuc_label, cell_label, single_nuc=False, cell_alone=False)
        path = os.path.join(path_output, filename+"_masks")
        plot.plot_images([nuc_label, cell_label], titles=["Labelled nuclei", "Labelled cells"], path_output=path)

        #decompose dense regions
        if run_all or input_thresholds[manual_counter, 1] != reference_threshold[manual_counter, 1]:
            RNA_1_spot_decomposition, dense_regions_1, reference_spot_1 = detection.decompose_dense(
            image=RNA_1_im, 
            spots=RNA_1_spot, 
            voxel_size=(300, 107.5, 107.5), 
            spot_radius=RNA_spot_radius, 
            alpha=alpha,  # alpha impacts the number of spots per candidate region
            beta=beta,  # beta impacts the number of candidate regions to decompose
            gamma=gamma)  # gamma the filtering step to denoise the image
            RNA_1_spot_IN, RNA_1_spot_Out = multistack.identify_objects_in_region(cell_label - nuc_label, RNA_1_spot, 3)
            RNA_1_spot_declustered_IN, RNA_1_spot_declustered_Out = multistack.identify_objects_in_region(cell_label - nuc_label, RNA_1_spot_decomposition, 3)
            path_RNA_1 = os.path.join(path_output, filename + '_mRNA_1_spots.csv')
            path_RNA_1_ind = os.path.join(path_output, filename + '_mRNA_1_spots_ind.csv')
            stack.save_data_to_csv(RNA_1_spot_declustered_IN, path_RNA_1, delimiter=',') 
            stack.save_data_to_csv(RNA_1_spot_IN, path_RNA_1_ind, delimiter=',')
            path_1 = os.path.join(path_output, RNA_1_filename+"_declustered_spots")
            #plot mRNA spots after decomposition of dense regions
            plot.plot_detection(RNA_1_mip, RNA_1_spot_declustered_IN, title=RNA_1_filename+" Decomposed dense regions", contrast=True, path_output=path_1)
            plot.plot_detection(RNA_1_mip, RNA_1_spot_IN, title=RNA_1_filename, contrast=True)

        # if run_all or input_thresholds[manual_counter, 2] != reference_threshold[manual_counter, 2]:
        #     RNA_2_spot_decomposition, dense_regions_2, reference_spot_2 = detection.decompose_dense(
        #     image=RNA_2_im, 
        #     spots=RNA_2_spot, 
        #     voxel_size=(300, 107.5, 107.5), 
        #     spot_radius=RNA_spot_radius, 
        #     alpha=alpha,  # alpha impacts the number of spots per candidate region
        #     beta=beta,  # beta impacts the number of candidate regions to decompose
        #     gamma=gamma)  # gamma the filtering step to denoise the image
        #     RNA_2_spot_IN, RNA_2_spot_Out = multistack.identify_objects_in_region(cell_label - nuc_label, RNA_2_spot, 3)
        #     RNA_2_spot_declustered_IN, RNA_2_spot_declustered_Out = multistack.identify_objects_in_region(cell_label - nuc_label, RNA_2_spot_decomposition, 3)
        #     path_RNA_2 = os.path.join(path_output, filename + '_mRNA_2_spots.csv')
        #     path_RNA_2_ind = os.path.join(path_output, filename + '_mRNA_2_spots_ind.csv')
        #     stack.save_data_to_csv(RNA_2_spot_declustered_IN, path_RNA_2, delimiter=',') 
        #     stack.save_data_to_csv(RNA_2_spot_IN, path_RNA_2_ind, delimiter=',') 
        #     path_2 = os.path.join(path_output, RNA_2_filename+"_declustered_spots")
        #     plot.plot_detection(RNA_2_mip, RNA_2_spot_declustered_IN, title=RNA_2_filename+" Decomposed dense regions", contrast=True, path_output=path_2)
        #     plot.plot_detection(RNA_2_mip, RNA_2_spot_IN, title=RNA_2_filename, contrast=True)

        if run_all or input_thresholds[manual_counter, 0] != reference_threshold[manual_counter, 0]:
            IF_spot_IN, IF_spot_OUT = multistack.identify_objects_in_region(cell_label - nuc_label, IF_spot, 3)
            path_IF = os.path.join(path_output, filename + '_PB_spots.csv')
            stack.save_data_to_csv(IF_spot_IN, path_IF, delimiter=',')
            path = os.path.join(path_output, IFfilename+"_spots")
            plot.plot_detection(IF_mip, IF_spot_IN, radius=4, title=IFfilename+"Detected spots", contrast=True, path_output=path)
    

headers = ['File Name', 'P-Body Threshold', "RNA1 Threshold"]
result_file_name = "HSP90AA+HSP110 Rerun.csv"
if os.path.isfile(result_file_name):
    user_input = input("Results file already exists do you want to overwrite(y/n): ")
    if user_input == "y":
        write_csv(result_file_name, manual_thresholds, headers)
    else:
        user_input = input("What is the new results file name (must end in .csv): ")
        write_csv(user_input, manual_thresholds, headers)
else:
    write_csv(result_file_name, manual_thresholds, headers)