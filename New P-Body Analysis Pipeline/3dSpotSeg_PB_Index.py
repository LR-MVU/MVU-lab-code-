
"""
Code to extract size data from P-Body masks and colocalization data between P-Body and 
mRNA masks.
 
@author: Lydia Hodgins
Affilitation: Vera Ugalde Lab, McGill University
"""
#################
#Import Packages
#################
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import pandas as pd
import os
from bigfishUI import *
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.classification as classification
from tifffile import imread
from skimage import measure
from skimage.measure import regionprops
from skimage.draw import disk
import math
from scipy import ndimage
from imageio import imread as imread2
from scipy.ndimage import label
from scipy.ndimage import binary_dilation
import matplotlib.pyplot as plt

##########
#Functions
##########
def get_filenames(directory, chan, filetype = '.tif'):
    '''gets all filenames corresponding to input channel
    str, str -> 2D list of match objects'''
    
    lis = os.listdir(directory)
    
    tif_pattern = re.compile(".*" + chan + ".*" + filetype)
    tif_files = [f for f in os.listdir(directory) if re.match(tif_pattern, f)]

    tif_files.sort()
    #print(tif_files)
    return tif_files 

def localized_PB_Index(PB_mask, rna_mask, PBimg, rnaimg, cell_area, pixel_dialation = 2):
    ''' computes the P-Body index for each P-Body considering only an area 2r^2 around the P-Body
    returns the list of individual P-Body indicies and average cell index 
    '''
    
    PBs = np.unique(PB_mask) #get list of individual PBody indices
   
    Int_INs = [] #list for intensity within the P-Body
    Int_OUTs = [] #list for intensity outside of the P-Body
    par_coef = np.zeros(len(PBs)-1) # Array to hold partition coefficients
    PB_Indices = np.zeros(len(PBs)-1) # Array to hold P-Body Indices
    PB_fraction = np.zeros(len(PBs)-1) # Array to hold P-Body Indices
    count = 0

    #Loop through all P-Bodies
    for j in PBs[1:]:
        zs = []
        zmax = 0
        Amax = 0

        #Find focal plane (zmax)
        for z in range(PB_mask.shape[0]):
            if j in PB_mask[z,:,:]:
                zs.append(z)
                area = np.count_nonzero(PB_mask[z,:,:] == j)
                if area > Amax: 
                    zmax = z
                    Amax = area
        
        #Dialte the P-Body mask to ensure OUT concentration is not too close to PB interface
        temp_mask = np.zeros((PB_mask.shape[1], PB_mask.shape[2]))
        temp_mask[PB_mask[zmax,:,:]==j] = 1
        structure = np.ones((3, 3), dtype=np.uint8)
        dialate_mask = ndimage.binary_dilation(temp_mask, structure=structure).astype(temp_mask.dtype)
        dialate_mask = ndimage.binary_dilation(dialate_mask, structure=structure).astype(dialate_mask.dtype)
        
        #Set PB center to position with highest pixel intensity
        PB = PBimg[zmax,:,:]*dialate_mask
        x = np.unravel_index(np.argmax(PB), PB.shape)[0]
        y = np.unravel_index(np.argmax(PB), PB.shape)[1]

        #Compute the radius of P-Body and make mask for ring around P-Body
        props1 = regionprops(dialate_mask.astype(int))
        if props1:
            r= int(2*math.ceil(props1[0].axis_major_length))
        else: r = 0
        if r != 0:
            local_PB = PB_mask[zmax-1:zmax+2, x-r:x+r+1, y-r:y+r+1]# P-Body mask in local area
            dialate = dialate_mask[x-r:x+r+1, y-r:y+r+1]
            dialate = np.repeat(dialate[np.newaxis,...], 3, axis=0) 
            local_mRNA = rna_mask[zmax-1:zmax+2, x-r:x+r+1, y-r:y+r+1] #RNA mask in local area
            rna_crop = rnaimg[zmax-1:zmax+2,x-r:x+r+1, y-r:y+r+1] #RNA image in local area

            #Check that P-Body is not too close to any edge of the image
            if (np.shape(local_PB)[0] == 3) and (np.shape(local_PB)[1] == 2*r+1) and (np.shape(local_PB)[2] == 2*r+1):
                #Make mask for outside ring
                OUTmask = np.zeros(np.shape(local_PB))
                rr, cc = disk((r-1,r-1), r)
                OUTmask[:, rr, cc] = j
                OUTmask[(local_PB != 0)] = 0
                OUTmask[(dialate != 0)] = 0
                local_PB = (dialate.astype(int)) * j
                n_IN = sum(local_mRNA[local_PB==j]) #Number of mRNA inside PB
                n_OUT = sum(local_mRNA[OUTmask==j]) #Number of mRNA outside PB
                
                PB_vol = len(local_PB[local_PB==j]) #Volume of PB
                cyt_vol = len(local_PB[OUTmask==j]) #Volume of considered cytoplasm

                #Compute mean RNA intensity inside and outside of PB
                props = regionprops(local_PB, rna_crop)
                Int_IN = props[0].intensity_mean
                props = regionprops(OUTmask.astype(int), rna_crop)
                Int_OUT = props[0].intensity_mean
                Int_INs.append(Int_IN)
                Int_OUTs.append(Int_OUT)

                #Compute PB Index and PB partition coefficient
                if n_IN + n_OUT != 0:
                    index = n_IN / ((n_IN + n_OUT)*(PB_vol/(PB_vol + cyt_vol)))
                    PB_Indices[count] = (index)
                    coef = (n_IN/PB_vol) / (n_OUT/cyt_vol)
                    par_coef[count] = (coef)
                    fraction =  n_IN/(n_IN+n_OUT)
                    PB_fraction[count] = (fraction)
                else:
                    PB_Indices[count] = (np.nan)
                    par_coef[count] = (np.nan)
                    PB_fraction[count] = (np.nan)
            else:
                PB_Indices[count] = (np.nan)
                par_coef[count] = (np.nan)
                PB_fraction[count] = (np.nan)
        else:
            PB_Indices[count] = (np.nan)
            par_coef[count] = (np.nan)
            PB_fraction[count] = (np.nan)
        count += 1
    avg_Index = np.nanmean(PB_Indices) #Average P-Body index for the cell
    
    # Adjusted part for P-Body Index
    binary_pb_array = np.where(PB_mask > 0, 1, 0)
    dialated_pb_mask = np.where(binary_dilation(binary_pb_array, iterations = pixel_dialation) == True, 1, 0)
    mrna_in_pb_mask = dialated_pb_mask * rna_mask
    num_rna_in_pb = np.sum(mrna_in_pb_mask)
    total_rna = np.sum(rna_mask)
    twoD_binary_pb_mask = np.max(binary_pb_array, axis=2)
    total_pb_area = np.sum(twoD_binary_pb_mask)
    pb_index_total_cell = num_rna_in_pb / (total_rna * (total_pb_area / cell_area))
    
    return avg_Index, PB_Indices, par_coef, Int_INs, Int_OUTs, PB_fraction, pb_index_total_cell

def pbody_coords_optimized(pbody_mask_array):
    pbody_coord_list = []
    image = np.zeros((pbody_mask_array.shape[0],pbody_mask_array.shape[1],3), dtype = np.uint8)
    binary_mask = pbody_mask_array > 0
    array_labeled, num_pbodies = label(binary_mask)
    
    for i in range(1, num_pbodies + 1):
        coords = np.argwhere(array_labeled == i)
        pbody_coord_list.append(tuple(coords[0]))
    """    
    for i in range(len(pbody_coord_list)):
        y = pbody_coord_list[i][0]
        x = pbody_coord_list[i][1]
        image[max(0,y-4):min(pbody_mask_array.shape[0],y+4), max(0,x-4):min(pbody_mask_array.shape[1],x+4) ] = (255,0,0)
        
    plt.imshow(image)
    plt.axis('off')
    plt.show()
    """
    
    return np.array(pbody_coord_list)
        
        
        
                


##############################
# Start of Analysis
###############################

#Get root directory from user
path_root = getImagesDir() 

Conditions = [
#   "ActD+HS_3HR", 
#   "ActD+HS_6HR", 
#   "ActD+HS_9HR", 
#   "ActD+HS_12HR", 
# "ActD+HS_HS",
#   "CT_CT",  
#   "CT_3HR",
#   "CT_6HR",
#   "CT_9HR", 
#   "CT_12HR",
#   "HS+ActD_3HR",
#   "HS+ActD_6HR",
#   "HS+ActD_9HR",
#   "HS+ActD_12HR",
#   "HS+ActD_HS"
                   'Ctl', 
                     'HS', 
                     'HS_2H_R', 
                   'HS_3H_R', 
                      'HS_4H_R', 
                       'HS_6H_R', 
                     'HS_8H_R', 
                    'HS_10H_R', 
                     'HS_12H_R',
              #  'HS_24H_R'
             #'12H_Ctl_ActD_Rep2',
             #'12H_Ctl_ActD_rep3'
             #'MEFs_Dcp1a'
               # 'WT_CT',
               # 'WT_HS', 
               # 'WT_2HR', 
               # 'WT_3HR', 
               # 'WT_4HR', 
               # 'WT_6HR', 
               # 'WT_8HR', 
               # 'HNR_CT', 
               # 'HNR_HS',
               # 'HNR_2HR',
               # 'HNR_3HR',
               # 'HNR_4HR',
               # 'HNR_6HR',
               # 'HNR_8HR',
               #  'CT',
                # '24HR'
              ]
dx = 0.1075 #um (pixel resolution)
PBsize_dict = {} #dictionary to hold all PB sizes
PBInd_dict = {} #dictionary to hold all PB Indices
wholePBInd_dict = {}
PB_R_dict = {} #dictionary to hold all PB radii
PBInt_dict = {} #dictionary to hold all PB integrated intensities
PBmaxInt_dict = {} #dictionary to hold all PB max intensities
numPB_dict = {} #dictionary to hold # of PBs/cell
RNAInt_dict = {} #dictionary to hold RNA intensities
Inti_dict = {} #dictionary to hold RNA intensity inside PB
Into_dict = {} #dictionary to hold RNA intensity outside PB
par_coef_dict = {} #dictionary to hold PB partition coefficient
RNAnum_dict = {} #dictionary to hold # RNA/cell
PB_fraction_dict = {}
avg_loc_index_dict = {}

#Loop through all conditions
for C in Conditions:
    print(C)
    path_input = os.path.join(path_root, C)
    path_PB = os.path.join(path_input, '3dSpotSeg') #PB masks stored in 3dSpotSeg/
    #get lists of files
    cell_mask_file = get_filenames(path_input, 'Print')
    nuc_mask_file = get_filenames(path_input, 'nuc')
    PB_mask_file = get_filenames(path_PB, '640S')
    PB_im_file = get_filenames(path_input, '640S')
    max_PB_mask_file = get_filenames(path_PB, 'ax')
    rna_im_file = get_filenames(path_input, '555S')

    PB_mask_file = list(filter(lambda x: 'max' not in x, PB_mask_file))
    
    #Inialize lists to store data for the condition
    PB_areas = [] # pxls
    PB_Rs = [] # um 
    PB_Int = [] #integrated intensity
    PB_max_int = [] #max intensity
    PB_Inds = []
    whole_PB_inds = []
    PB_num = []
    Int_Inds = []
    RNA_Int = [] #mRNA integrate intensity
    Int_Is = []
    Int_Os = []
    par_coefs = []
    RNAnum = []
    PB_fraction = []
    PB_Inds_cell = []

    #Loop through each cell
    for i in range(len(PB_mask_file)):
        print(PB_mask_file[i])
        #Load data
        PB_mask = imread(os.path.join(path_PB, PB_mask_file[i]))
        PB_mask_mask = imread(os.path.join(path_PB, max_PB_mask_file[i]))
        PBimg = imread(os.path.join(path_input, PB_im_file[i]))
        RNAimg = imread(os.path.join(path_input, rna_im_file[i]))
        smFISH_mip = stack.maximum_projection(RNAimg)
        cell_mask  = plt.imread(os.path.join(path_input, cell_mask_file[i]))
        nuc_mask = plt.imread(os.path.join(path_input, nuc_mask_file[i]))
        rna_file = cell_mask_file[i]
        rna_file = rna_file[:-4]
        try:
            RNA_locs = pd.read_csv(os.path.join(path_input, 'output', rna_file + '_mRNA_spots.csv'), sep=',', header=None)
            if RNA_locs.empty:
                print(rna_file, 'is empty')
            else:
                RNA_locs=RNA_locs.values
                print(rna_file + '_mRNA_1_spots.csv')
                
                #create mRNA mask
                RNA_mask = np.zeros(np.shape(PBimg))
                for j in range(RNA_locs.shape[0]):
                    z = RNA_locs[j][0]
                    y = RNA_locs[j][1]
                    x = RNA_locs[j][2]
                    RNA_mask[z,y,x] += 1
            
        except pd.errors.EmptyDataError:           
            print(rna_file, 'is empty')

            RNA_mask = np.zeros(np.shape(PB_mask))

        #Remove PBs and RNA that are outside cell or in nucleus
        cyt_mask = cell_mask - nuc_mask
        PB_mask_mask = PB_mask_mask * cyt_mask
        RNA_mask = RNA_mask * cyt_mask
        RNAnum.append(np.sum(RNA_mask))
        PBs = np.unique(PB_mask_mask)
        PB_num.append(len(PBs[1:]))
         
        #Loop through each PB
        for ID in PBs[1:]:
            #Compute radius and Integrated Intensity of PB
            area = []

            for z in range(PB_mask.shape[0]):
                area.append(np.count_nonzero(PB_mask[z,:,:] == ID))
            PB_areas.append(max(area))
            zmax = np.argmax(np.array(area))
            
            
            #print(pb_coord)
            if zmax > 0 and zmax < 39:
                slice = measure.label(PB_mask[zmax-1:zmax+2,:,:] == ID)
                props = regionprops(slice, PBimg[zmax-1:zmax+2,:, :])
                PB_Rs.append((props[0].equivalent_diameter_area)*dx/2)
                PB_Int.append((props[0].area)*(props[0].intensity_mean))
                PB_max_int.append((props[0].intensity_max))
                props2 = regionprops(slice, RNAimg[zmax-1:zmax+2,:,:])
                RNA_Int.append((props2[0].area)*(props2[0].intensity_mean))
            else:
                slice = measure.label(PB_mask[zmax,:,:] == ID)
                props = regionprops(slice, PBimg[zmax,:,:])
                PB_Rs.append((props[0].equivalent_diameter_area)*dx/2)
                PB_Int.append(props[0].intensity_mean)
                PB_max_int.append(props[0].intensity_max)
                props2 = regionprops(slice, RNAimg[zmax,:,:])
                RNA_Int.append((props2[0].area)*(props2[0].intensity_mean))

        
        #Get colocalization measurements
        
        RNA_locs = np.int32(RNA_locs)
        pb_coords = pbody_coords_optimized(PB_mask_mask)
        
        #### PB-index functions from BigFish #####
        input_features = classification.prepare_extracted_data(cell_mask=cell_mask-nuc_mask,
                                                          nuc_mask=nuc_mask,
                                                          ndim = 3,
                                                          rna_coord=RNA_locs,
                                                          centrosome_coord = pb_coords)
        
        # Returns cell area without nucleus
        area = classification.features_area(input_features[0], input_features[5], input_features[6])[3]
        
       
        avg_loc_Index, Indices, par_coef, Int_I, Int_O, PB_fraction, total_pbody_index = localized_PB_Index(PB_mask, RNA_mask, PBimg, RNAimg, area)
        PB_Inds += Indices.tolist()
        Int_Is += Int_I
        Int_Os += Int_O
        par_coefs += par_coef.tolist()
        PB_fraction += PB_fraction.tolist()
        PB_Inds_cell.append(avg_loc_Index)
        whole_PB_inds.append(total_pbody_index)

    #Add data to dictionary with condition as key
    PBsize_dict[C] = PB_areas
    PBInd_dict[C] = PB_Inds
    PB_R_dict[C] = PB_Rs
    PBInt_dict[C] = PB_Int
    PBmaxInt_dict[C] = PB_max_int
    numPB_dict[C] = PB_num
    RNAInt_dict[C] = RNA_Int
    Inti_dict[C] = Int_Is
    Into_dict[C] = Int_Os
    par_coef_dict[C] = par_coefs
    RNAnum_dict[C] = RNAnum
    PB_fraction_dict[C] = PB_fraction
    avg_loc_index_dict[C] = PB_Inds_cell
    wholePBInd_dict[C] = whole_PB_inds
    
#Save all dictionaries to .csv files
# size_df = pd.DataFrame.from_dict(PBsize_dict, orient='index')
# size_df = size_df.transpose()
# size_df.to_csv(os.path.join(path_root, 'PBsize_v_cond.csv'))

# Ind_df = pd.DataFrame.from_dict(PBInd_dict, orient='index')
# Ind_df = Ind_df.transpose()
# Ind_df.to_csv(os.path.join(path_root, 'PBIndex_v_cond.csv')) 

# Ind_cell_df = pd.DataFrame.from_dict(avg_loc_index_dict, orient='index')
# Ind_cell_df = Ind_cell_df.transpose()
# Ind_cell_df.to_csv(os.path.join(path_root, 'PBIndex_cell_v_cond.csv')) 

wholeInd_cell_df = pd.DataFrame.from_dict(wholePBInd_dict, orient='index')
wholeInd_cell_df = wholeInd_cell_df.transpose()
wholeInd_cell_df.to_csv(os.path.join(path_root, 'wholePBIndex_cell_v_cond.csv')) 

# R_df = pd.DataFrame.from_dict(PB_R_dict, orient='index')
# R_df = R_df.transpose()
# R_df.to_csv(os.path.join(path_root, 'PBradius_v_cond.csv'))

# Int_df = pd.DataFrame.from_dict(PBInt_dict, orient='index')
# Int_df = Int_df.transpose()
# Int_df.to_csv(os.path.join(path_root, 'PBIntmean_v_cond.csv'))

# maxInt_df = pd.DataFrame.from_dict(PBmaxInt_dict, orient='index')
# maxInt_df = maxInt_df.transpose()
# maxInt_df.to_csv(os.path.join(path_root, 'PBmaxInt_v_cond.csv'))

# numPB_df = pd.DataFrame.from_dict(numPB_dict, orient='index')
# numPB_df = numPB_df.transpose()
# numPB_df.to_csv(os.path.join(path_root, 'numPBs_v_cond.csv'))

# numRNA_df = pd.DataFrame.from_dict(RNAnum_dict, orient='index')
# numRNA_df = numRNA_df.transpose()
# numRNA_df.to_csv(os.path.join(path_root, 'numRNA_v_cond.csv'))

# RNAInt_df = pd.DataFrame.from_dict(RNAInt_dict, orient='index')
# RNAInt_df = RNAInt_df.transpose()
# RNAInt_df.to_csv(os.path.join(path_root, 'RNAIntensity_v_cond.csv'))

# IntIn_df = pd.DataFrame.from_dict(Inti_dict, orient='index')
# IntIn_df = IntIn_df.transpose()
# IntIn_df.to_csv(os.path.join(path_root, 'IntensityIn_v_cond.csv'))

# IntOut_df = pd.DataFrame.from_dict(Into_dict, orient='index')
# IntOut_df = IntOut_df.transpose()
# IntOut_df.to_csv(os.path.join(path_root, 'IntensityOut_v_cond.csv'))

# par_coef_df = pd.DataFrame.from_dict(par_coef_dict, orient='index')
# par_coef_df = par_coef_df.transpose()
# par_coef_df.to_csv(os.path.join(path_root, 'ParCoef_v_cond.csv'))

# PB_fraction_df = pd.DataFrame.from_dict(PB_fraction_dict, orient='index')
# PB_fraction_df = PB_fraction_df.transpose()
# PB_fraction_df.to_csv(os.path.join(path_root, 'PBfraction_v_cond.csv'))