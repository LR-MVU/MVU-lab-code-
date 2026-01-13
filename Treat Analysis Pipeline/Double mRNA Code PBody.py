import os
import sys
from skimage import io
import numpy as np
import pandas as pd
import csv
from natsort import natsorted
import matplotlib.pyplot as plt
import math

# The write_csv function takes a file_location as the file_name, an array of data, and a list of headers
# and produces a csv file
def write_csv(file_name, array, headers):
    with open(file_name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(array)

def csv_reader (csv_folder, csv_file):
    try:
        df_csv = pd.read_csv(os.path.join(csv_folder, csv_file), dtype = str)
        return df_csv.to_numpy()
    except pd.errors.EmptyDataError:
        return np.array([])

efficiency_chan1 = 0.9
efficiency_chan2 = 0.9
folder_channel1 = r"C:\Users\ugald\Downloads\double mRNA\u2os transfection output bigFISH\mRNA_1"
folder_channel2 = r"C:\Users\ugald\Downloads\double mRNA\u2os transfection output bigFISH\mRNA_2"
folder_pbody = r""
chan1 = "mRNA_1"
chan2 = "mRNA_2"
pbody = ""
pixel_size = 107.5
chromo_abberation = 3

channel1_np_array = {}
for file in os.listdir(folder_channel1):
    if os.path.isfile(os.path.join(folder_channel1, file)) and "csv" in file and "spots" in file:
        channel1_np_array[file.split(chan1)[0]] = csv_reader(folder_channel1, file)

channel2_np_array = {}
for file in os.listdir(folder_channel2):
    if os.path.isfile(os.path.join(folder_channel2, file)) and "csv" in file and "spots" in file:
        channel2_np_array[file.split(chan2)[0]] = csv_reader(folder_channel2, file)

pbody_np_array = {}
for file in os.listdir(folder_pbody):
    if os.path.isfile(os.path.join(folder_pbody, file)) and "csv" in file and "spots" in file:
        pbody_np_array[file.split(pbody)[0]] = csv_reader(folder_pbody, file)

def distance(coord1, coord2):
    return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

def extract_coords(sheet):
    coord_list = []
    for i in range(sheet.shape[0]):
        z_coord = int(sheet[i,0])
        y_coord = int(sheet[i,1])
        x_coord = int(sheet[i,2])
        coord_list.append((z_coord, y_coord, x_coord)) # (y,x)
    return coord_list

def get_max_z(mRNA_1_coords, mRNA_2_coords, pbody_coords):
    max_z = 0 
    for i in mRNA_1_coords:
        if i[0] > max_z:
            max_z = i[0]
    for i in mRNA_2_coords:
        if i[0] > max_z:
            max_z = i[0]
    for i in pbody_coords:
        if i[0] > max_z:
            max_z = i[0]
    return max_z

def create_pbody_mask(coords, max_z, pbody_dialation = 2):
    pbody_mask = np.zeros((2304, 2304, max_z), dtype = int)
    for i in range(len(coords)):
        x, y, z = coords[i][2]-1, coords[i][1]-1, coords[i][0]-1
        pbody_mask[max((x-pbody_dialation),0):min((x+pbody_dialation+1),2303), max((y-pbody_dialation),0):min((y+pbody_dialation)+1,2303), max((z-pbody_dialation),0):min((z+pbody_dialation+1),2303)] = 1
    return pbody_mask

def mRNA_cluster(mRNA_coords, max_z, cube_size = 2):
    cluster_image = np.zeros((2304, 2304, max_z), dtype = int)
    for i in range(len(mRNA_coords)):
        x, y, z = mRNA_coords[i][2]-1, mRNA_coords[i][1]-1, mRNA_coords[i][0]-1
        for j in range(len(mRNA_coords)):
            if distance(mRNA_coords[i], mRNA_coords[j]) < chromo_abberation and i != j:
                cluster_image[max((x-cube_size),0):min((x+cube_size+1),2303), max((y-cube_size),0):min((y+cube_size)+1,2303), max((z-cube_size),0):min((z+cube_size+1),2303)] = 1
    return cluster_image

def split_mRNA(cluster_mask, coords):
    disperse_mrna = []
    clustered_mrna = []
    for coord in coords:
        x = coord[2] - 1
        y = coord[1] - 1
        z = coord[0] - 1
        if cluster_mask[x,y,z] == 1:
            clustered_mrna.append(coord)
        else:
            disperse_mrna.append(coord)
    return disperse_mrna, clustered_mrna    

            
def pair_coords(chan1_coords, chan2_coords):
    pair_list = []
    distance_all = []
    distance_paired = []
    unpaired_chan1 = []
    for coord1 in chan1_coords:
        min_dist = sys.maxsize
        current_coord = (-1,-1,-1)
        for coord2 in chan2_coords:
            dist = distance(coord1, coord2)
            if min_dist > dist:
                min_dist = dist
                current_coord = coord2
        if min_dist <= chromo_abberation:
            pair_list.append((coord1,current_coord))
            chan2_coords.remove(current_coord)
            distance_paired.append(min_dist)
        else:
            unpaired_chan1.append(coord1)
        distance_all.append(min_dist)
    return pair_list, distance_all, distance_paired, unpaired_chan1, chan2_coords

def pbody_loc(mRNA_coords, pbody_mask):
    mRNA_pbody = []
    mRNA_no_pbody = []
    for coord in mRNA_coords:
        x = coord[2] - 1
        y = coord[1] - 1
        z = coord[0] - 1
        if pbody_mask[x,y,z] > 0:
            mRNA_pbody.append(coord)
        else:
            mRNA_no_pbody.append(coord)
    return len(mRNA_pbody)/len(mRNA_coords)

def pbody_analysis(mRNA_list, pbody_mask):
    result = []
    for i in range(len(mRNA_list)):
        result.append(pbody_loc(mRNA_list[i], pbody_mask))
    return result
    
                
# Remove circles of area, where conflicts are on each channel and then use remaining mRNA to get a fraction paired
# In areas of clumping look at the total frequency or total amplitude of each channel
# Fraction of mRNA discarded vs the total number of mRNAs
# Look in the Z Stack, 200nm per Z

def visualization(imgname, cluster_mask, c1, c2):
    cluster_mask = np.argmax(cluster_mask, axis = 2).T
    image = np.dstack([cluster_mask*0]*3).astype(np.uint8)
    for _, y ,x in c1:
        image[y, x, 0] = 255 
    for _, y, x in c2:
        image[y, x, 2] = 255
    title = imgname
    plt.figure(figsize=(10, 10))
    plt.title(title)
    plt.imshow(image)
    plt.axis('off')
    plt.savefig(imgname, format = "tiff", dpi = 300)
    plt.show()

def histogram(data, title):
    plt.hist(data, bins=30, edgecolor='black')
    plt.xlabel("MS2-SunTag Distance")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.show()

final_results = []
total_distance_all = []
total_distance_paired = []
for key in list(channel1_np_array.keys()):
    if key in channel2_np_array and key in pbody_np_array:
        chan1_data = channel1_np_array.pop(key)
        chan2_data = channel2_np_array.pop(key)
        pbody_data = pbody_np_array.pop(key)
        chan1_coords = extract_coords(chan1_data)
        chan2_coords = extract_coords(chan2_data)
        pbody_coords = extract_coords(pbody_data)
        print("extracted coords")
        max_z = get_max_z(chan1_coords, chan2_coords, pbody_coords)
        cluster_array1 = mRNA_cluster(chan1_coords, max_z)
        cluster_array2 = mRNA_cluster(chan2_coords, max_z)
        pbody_mask = create_pbody_mask(pbody_coords, max_z)
        print("cluster masks created")
        merged_cluster_array = np.maximum(cluster_array1, cluster_array2)
        chan1_disperse_mrna, chan1_clustered_mrna = split_mRNA(merged_cluster_array, chan1_coords)
        chan2_disperse_mrna, chan2_clustered_mrna = split_mRNA(merged_cluster_array, chan2_coords)
        print("mRNA split")
        copy_chan2_disperse_mrna = chan2_disperse_mrna.copy()
        paired_coords, distance_all, distance_paired, unpaired_chan1, unpaired_chan2 = pair_coords(chan1_disperse_mrna, copy_chan2_disperse_mrna)
        pbody_results = pbody_analysis([unpaired_chan1, unpaired_chan2, paired_coords, chan1_clustered_mrna, chan2_clustered_mrna], pbody_mask)
        """
        total_distance_all.extend(distance_all)
        total_distance_paired.extend(distance_paired)
        print("mRNA paired")
        visualization(key+"clustered", merged_cluster_array, chan1_clustered_mrna, chan2_clustered_mrna)
        visualization(key+"disperse", merged_cluster_array,chan1_disperse_mrna, chan2_disperse_mrna)
        cluster_mask = np.argmax(merged_cluster_array, axis = 2).T
        image = np.dstack([cluster_mask*255]*3).astype(np.uint8)
        plt.figure(figsize=(10, 10))
        plt.imshow(image)
        plt.axis('off')
        plt.savefig("Background"+key, format = "tiff", dpi = 300)
        plt.show()
        """
        final_results.append([key, len(paired_coords), len(chan1_disperse_mrna)-len(paired_coords), len(chan2_disperse_mrna)-len(paired_coords), len(chan1_clustered_mrna), len(chan2_clustered_mrna)] + pbody_results)

#histogram(total_distance_all,"Distance Between MS2 and SunTag in Disperse mRNAs")
#histogram(total_distance_paired,"Distance Between MS2 and SunTag in Paired Disperse mRNAs")

headers = ["Filename", "Paired Disperse Signals", "Chan1 Unpaired Disperse Signal", "Chan2 Unpaired Disperse Signal", "Chan1 Frequency Clustered SIgnal", "Chan2 Freuqency Clustered Signal", "Unpaired Chan1 Pbody Loc", "Unpaired Chan2 Pbody Loc", "Paired Pbody Loc", "Clustered 1 Pbody Loc", "Clustered 2 Pbody Loc"]
file_name = "Results.csv"
write_csv(file_name, final_results, headers)



