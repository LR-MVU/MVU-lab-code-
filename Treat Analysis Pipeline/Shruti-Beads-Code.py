import os
import sys
from skimage import io
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import math

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

folder_channel1 = r"Z:\Charlotte\Codes From James\Total Beads\Spot detection results 555S automatic"
folder_channel2 = r"Z:\Charlotte\Codes From James\Total Beads\Spot detection results 640S automatic"
folder_prints = r"Z:\Charlotte\Codes From James\Total Beads\Total Skeleton"
chan1 = "555S"
chan2 = "640S"
pixel_size = 107.5

channel1_np_array = {}
for file in os.listdir(folder_channel1):
    if os.path.isfile(os.path.join(folder_channel1, file)) and "csv" in file and "bead" in file:
        channel1_np_array[file.split(chan1)[0]] = csv_reader(folder_channel1, file)

channel2_np_array = {}
for file in os.listdir(folder_channel2):
    if os.path.isfile(os.path.join(folder_channel2, file)) and "csv" in file and "bead" in file:
        channel2_np_array[file.split(chan2)[0]] = csv_reader(folder_channel2, file)
        
prints_np_array = []
for file in os.listdir(folder_prints):
    if os.path.isfile(os.path.join(folder_prints, file)) and "print" in file:
        image = io.imread(os.path.join(folder_prints,file)).astype(np.uint8)
        grey_scale_image = np.mean(image, axis=-1)
        final_image = (np.where(grey_scale_image > 172, 255, 0)).squeeze() # Yellow is 173 and Purple is 51
        prints_np_array.append(final_image)

def distance(coord1, coord2):
    return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)

def distance_2d(coord1, coord2):
    if len(coord1) > 2:
        coord1 = coord1[1:]
    if len(coord2) > 2:
        coord2 = coord2[1:]
    return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2)

def max_distance(image):
    shape = image.shape
    median_coords = (shape[0]/2, shape[1]/2) # (y,x)
    cell_pixels = np.argwhere(image > 0)
    max_distacne_from_middle = 0
    for pixel in cell_pixels:
        distance_value = distance_2d(median_coords, pixel)
        if distance_value > max_distacne_from_middle:
            max_distacne_from_middle = distance_value
    return max_distacne_from_middle, median_coords

def max_distance_list(image_list):
    max_distance_from_middle = 0
    median_coords_list = []
    for image in image_list:
        distance, median_coords = max_distance(image)
        median_coords_list.append(median_coords)
        if distance > max_distance_from_middle:
            max_distance_from_middle = distance
    if median_coords_list.count(median_coords_list[0]) == len(median_coords_list):
        return max_distance_from_middle, median_coords_list[0]
    else:
        raise "Images are of different sizes, middle coordinates are different"

max_distance_from_middle, median_coords = max_distance_list(prints_np_array)
print("Max Distance from Center: ", max_distance_from_middle)

def extract_coords(sheet):
    coord_list = []
    for i in range(sheet.shape[0]):
        coord_array = sheet[i,0].split(";")
        z_coord = int(coord_array[0])
        y_coord = int(coord_array[1])
        x_coord = int(coord_array[2])
        coord_list.append((z_coord, y_coord, x_coord)) # (z,y,x)
    return coord_list

def pair_coords(chan1_coords, chan2_coords):
    paired_coords = {}
    for coord1 in chan1_coords:
        min_distance = sys.maxsize
        c = (-1,-1,-1)
        for coord2 in chan2_coords:
            distance_value = distance(coord1, coord2)
            if distance_value < min_distance:
                min_distance = distance_value
                c = coord2
        paired_coords[(coord1,c)] = min_distance
    return paired_coords


def filter_coords(max_distance_from_middle, median_coords, paired_coords_dict):
    included_chromatic_abberation_list = []
    for key in paired_coords_dict:
        coord1, coord2 = key
        distance_coord1 = distance_2d(median_coords, coord1)
        distance_coord2 = distance_2d(median_coords, coord2)
        if distance_coord1 <= max_distance_from_middle and distance_coord2 <= max_distance_from_middle:
            included_chromatic_abberation_list.append(paired_coords_dict[key])
    return included_chromatic_abberation_list

def histogram(data):
    plt.hist(data, bins=30, edgecolor='black')
    plt.xlabel("Chromatic Abberation in Nanometers")
    plt.ylabel("Frequency")
    plt.title("Chromatic Abberation")
    plt.show()

def chromatic_abberation(chan1_dict, chan2_dict):
    total_list_chromo_abberation_distances = []
    for key in list(chan1_dict.keys()):
        if key in chan2_dict:
            chan1_data = chan1_dict.pop(key)
            chan2_data = chan2_dict.pop(key)
            chan1_coords = extract_coords(chan1_data)
            chan2_coords = extract_coords(chan2_data)
            print("Coords Extracted")
            paired_coords = pair_coords(chan1_coords, chan2_coords)
            print("Coords Paired")
            filter_chromatic_abberation_distances = filter_coords(max_distance_from_middle, median_coords, paired_coords)
            total_list_chromo_abberation_distances.extend(filter_chromatic_abberation_distances)
    
    length_list = len(total_list_chromo_abberation_distances)
    sorted_chromo_abberation = (np.array(sorted(total_list_chromo_abberation_distances)))*pixel_size
    histogram(sorted_chromo_abberation[0:math.floor(length_list*90/100)])
    ninety_fith_percentile = sorted_chromo_abberation[math.floor(length_list*95/100)]
    ninety_percentile = sorted_chromo_abberation[math.floor(length_list*90/100)]
    seventy_fith_percentile = sorted_chromo_abberation[math.floor(length_list*75/100)]
    fifty_percentile = sorted_chromo_abberation[math.floor(length_list*50/100)]
    print("95th Percentile Chromatic Abberation: ", ninety_fith_percentile)
    print("90th Percentile Chromatic Abberation: ", ninety_percentile)
    print("75th Percentile Chromatic Abberation: ", seventy_fith_percentile)
    print("50th Percentile Chromatic Abberation: ", fifty_percentile)
    return ninety_fith_percentile

chromatic_abberation(channel1_np_array, channel2_np_array)

    










