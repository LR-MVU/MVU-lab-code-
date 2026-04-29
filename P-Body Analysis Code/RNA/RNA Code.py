# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:05:38 2026

@author: Lab
"""
import numpy as np
import os
import pandas as pd

import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.multistack as multistack
import bigfish.segmentation as segmentation
import bigfish.plot as plot
import re
import csv
from bigfishUI import *


# Function Creates a CSV file from an array of data and a list of headers
def write_csv(file_name, array, headers):
    with open(file_name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(array)

# Function reads in a CSV file and converts data to a numpy array        
def csv_reader (csv_path):
   df_csv = pd.read_csv((csv_path), dtype = str)
   return df_csv.to_numpy()

# Reads in files within a specific directory that contain a specific key word (channel)
def get_filenames(directory, chan):
    '''gets all filenames corresponding to input channel
    str, str -> 2D list of match objects'''
    
    lis = os.listdir(directory)
    
    tif_pattern = re.compile(".*_" + chan + ".*.tif")
    matches = regexMatchFilter(lis, tif_pattern)
    tif_files = [match.group() for match in matches]
    if chan == "555S":
        tif_files = [file for file in tif_files if not ("nuc" in file or "print" in file)]
    tif_files.sort()
    return tif_files 

# DECLUSTERING PARAMETERS
alpha=0.5  # alpha impacts the number of spots per candidate region
beta=1 # beta impacts the number of candidate regions to decompose
gamma=7  # gamma the filtering step to denoise the image
RNA_spot_radius = (350, 150, 150) # in nanometer (one value per dimension zyx)    
Voxel_size = (300, 107.5, 107.5) # in nanometer (one value per dimension zyx)

# Asks the user for the number of RNAs they wish to detect, and asks them to input key phrase (channel)
# that differentiates specific RNA files from other files in directory
num_rnas = int(input("Please enter the number of RNAs you wish to analyze: "))
if num_rnas == 1:
    rna1_chan = input("Input channel of RNA 1 (ex. 555S): ")
elif num_rnas == 2:
    rna1_chan = input("Input channel of RNA 1 (ex. 555S): ")
    rna2_chan = input("Input channel of RNA 2 (ex. 640S): ")
else:
    raise ValueError("Number of RNAs must be either 1 or 2")

# Input cell and mask channels so they can be properly read in
cell_mask_chan = input("Input cell mask channel (ex. print): ")
nuc_mask_chan = input("Input nuc mask channel (ex. nuc): ")

# Enter the folder containing the cropped images, including RNA files and cell and nucleus prints
folder = input(r"Input folder path containing cropped images: ")
# Input condition list, within the inputed folder above their should be subfolders for each of the conditions
conditions = input("Input list of conditions seperated by commas with no additional spaces: ").split(",")

# You must first perform manual detection, if you then want to adjust the thresholds for RNA detection because too many
# or too few spots are detected then adjust the thresholds in a copy of the produced spreadsheet and run manual detection
manual = bool(int(input("Are you doing automatic detection (input '0') or manual detection (input '1'): ")))
thresh_rna1 = {}
thresh_rna2 = {}
ref_rna1 = {}
ref_rna2 = {}
if manual:
    # Thresh path is the path to the excel file containing the thresholds you have altered
    thresh_path = input(r"Input path to CSV with manual thresholds: ")
    raw_thresh = csv_reader(thresh_path)
    # Ref path is the path to the escel file containing the original manual thresholds, these will be referenced to see what
    # analysis must be rerun
    ref_path = input(r"Input path to CSV wuth reference thresholds: ")
    raw_ref = csv_reader(ref_path)
    # Reads in altered thresholds and reference thresholds, only reruns spot detection for RNA thresholds that have been changed
    if num_rnas == 1:
        for i in range(raw_thresh.shape[0]):
            thresh_rna1[raw_thresh[i,0]] = float(raw_thresh[i,1])
        for i in range(raw_ref.shape[0]):
            ref_rna1[raw_ref[i,0]] = float(raw_ref[i,1])
    if num_rnas == 2:
        for i in range(raw_thresh.shape[0]):
            thresh_rna1[raw_thresh[i,0]] = float(raw_thresh[i,1])
            split_name = raw_thresh[i,0].split(rna1_chan)
            key_rna2 = split_name[0] + rna2_chan + split_name[1]
            thresh_rna2[key_rna2] = float(raw_thresh[i,2])
        for i in range(raw_ref.shape[0]):
            ref_rna1[raw_ref[i,0]] = float(raw_ref[i,1])
            split_name = raw_ref[i,0].split(rna1_chan)
            key_rna2 = split_name[0] + rna2_chan + split_name[1]
            ref_rna2[key_rna2] = float(raw_ref[i,2])
        
# Main function that performs RNA thresholding for a specific RNA channel
def rna_detection(RNA_chan, manual, thresh=None, ref=None):
    print("")
    print(RNA_chan)
    rna_thresh_list = []
    for cond in conditions:
        print(cond)
        path_input = os.path.join(folder, cond)
        path_output = os.path.join(path_input, "output")
        os.makedirs(path_output, exist_ok=True)
        
        # Gets lists of file names for RNA channel and nuclear and cell masks
        cell_mask_files = get_filenames(path_input, cell_mask_chan)
        nuc_mask_files = get_filenames(path_input, nuc_mask_chan)
        RNA_files = get_filenames(path_input, RNA_chan)
        
        for i in range(len(RNA_files)):
            RNA_filename = RNA_files[i][:-4]
            
            # Skip file if manual and the reference and inputed threshold are the same
            if manual:
                if RNA_filename in thresh and RNA_filename in ref:
                    if round(thresh[RNA_filename],3) == round(ref[RNA_filename], 3):
                        continue
            
            # Read in RNA image and create max projection
            TIF_path = os.path.join(path_input, RNA_files[i])
            RNA_im = stack.read_image(TIF_path)
            RNA_mip = stack.maximum_projection(RNA_im)
            
            # Produce an elbow plot of pixel intensity and save to output folder for QC
            plot.plot_elbow(
                images=RNA_im, 
                voxel_size=Voxel_size, 
                spot_radius=RNA_spot_radius,
                path_output=os.path.join(path_output, RNA_filename+"_elbow"))
            
            # Perform spot detection using automatic threshold and return and store threshold
            if not manual:
                RNA_spot, RNA_threshold = detection.detect_spots(images=RNA_im, return_threshold=True, 
                                                                     voxel_size=Voxel_size,  
                                                                     spot_radius=RNA_spot_radius)
                rna_thresh_list.append([RNA_filename, RNA_threshold])
                
            # Use manually inputed threshold to threshold RNA    
            if manual:
                if RNA_filename in thresh:
                    manual_thresh = thresh[RNA_filename]
                    
                    RNA_spot = detection.detect_spots(images=RNA_im, return_threshold=False, threshold = manual_thresh,
                                                                        voxel_size=Voxel_size,  
                                                                        spot_radius=RNA_spot_radius)
                else:
                    print(RNA_filename)
                    raise ValueError("You selected manual thresholding but one of your images did not have a threshold (above)")
                
            
            ###### Get Nucleus Mask #######
            TIF_path = os.path.join(path_input, nuc_mask_files[i])
            filename = nuc_mask_files[i][:-4]
            mask_im = stack.read_image(TIF_path)
            mask_im = mask_im.astype(bool)
            nuc_label = segmentation.label_instances(mask_im)
            
            ###### Get Cell Masks ######
            TIF_path = os.path.join(path_input, cell_mask_files[i])
            filename = cell_mask_files[i][:-4]
            mask_im = stack.read_image(TIF_path)
            mask_im = mask_im.astype(bool)
            cell_label = segmentation.label_instances(mask_im)
            
            # Match nucleus mask with cell mask
            nuc_label, cell_label = multistack.match_nuc_cell(nuc_label, cell_label, single_nuc=False, cell_alone=False)
            path = os.path.join(path_output, filename+"_masks")
            plot.plot_images([nuc_label, cell_label], titles=["Labelled nuclei", "Labelled cells"], path_output=path)
            
            
            #Perform RNA decomposition to esitmate the number of densly clustered RNAs
            RNA_spot_decomposition, dense_regions, reference_spot = detection.decompose_dense(
            image=RNA_im, 
            spots=RNA_spot, 
            voxel_size=Voxel_size, 
            spot_radius=RNA_spot_radius, 
            alpha=alpha,  # alpha impacts the number of spots per candidate region
            beta=beta,  # beta impacts the number of candidate regions to decompose
            gamma=gamma)  # gamma the filtering step to denoise the image
            RNA_spot_IN, RNA_spot_Out = multistack.identify_objects_in_region(cell_label - nuc_label, RNA_spot, 3)
            RNA_spot_declustered_IN, RNA_spot_declustered_Out = multistack.identify_objects_in_region(cell_label - nuc_label, RNA_spot_decomposition, 3)
            path_RNA = os.path.join(path_output, filename + '_mRNA' + RNA_chan + '_spots.csv')
            path_RNA_decomp = os.path.join(path_output, filename + '_mRNA' + RNA_chan + '_spots_decomp.csv')
            # Save both the decomposed and raw coordinates or detected RNAs in a excel spreadsheet
            stack.save_data_to_csv(RNA_spot_declustered_IN, path_RNA_decomp, delimiter=',') 
            stack.save_data_to_csv(RNA_spot_IN, path_RNA, delimiter=',')
            path_1 = os.path.join(path_output, RNA_filename+"_declustered_spots")
            #plot mRNA spots after decomposition of dense regions
            plot.plot_detection(RNA_mip, RNA_spot_declustered_IN, title=RNA_filename+" Decomposed dense regions", contrast=True, path_output=path_1)
            plot.plot_detection(RNA_mip, RNA_spot_IN, title=RNA_filename, contrast=True)
        
    return rna_thresh_list

# Executes the main function and calls write_csv to save the saved thresholds if automatic
if num_rnas == 1:
    rna1_thresh = rna_detection(rna1_chan, manual, thresh_rna1, ref_rna1)
    if not manual:
        headers = ["Filename", "RNA1 Threshold"]
        write_csv("RNA Thresholds.csv", rna1_thresh, headers)
    
elif num_rnas == 2:
    rna1_thresh = rna_detection(rna1_chan, manual, thresh_rna1, ref_rna1)
    rna2_thresh = rna_detection(rna2_chan, manual, thresh_rna2, ref_rna2)
    if not manual:
        rna1_thresh = np.array(rna1_thresh, dtype=object)
        rna2_thresh = np.array(rna2_thresh, dtype=object)
        headers = ["Filename", "RNA1 Threshold", "RNA2 Threshold"]
        rna_thresh = np.column_stack((rna1_thresh, rna2_thresh[:, 1]))
        write_csv("RNA Thresholds.csv", rna_thresh, headers)
            
        
        
        
        