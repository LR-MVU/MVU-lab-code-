#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 16:52:57 2023

@author: stellaxu
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 08:45:48 2023

@author: stellaxu
"""
# -*- coding: utf-8 -*-
"""
Created on Wed May 3 14:06:43 2023

@author: Lokha, Sethu
"""

import numpy as np
import math
import os
import re
import csv
import sys

# CSV Writer turns arrays into CSV files with given headers at specified location
def write_csv(filename, data, headers):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for row in data:
            writer.writerow(row)

# Determines the median amplitude of a single peptide by taking the median of all peptides not within an mRNA
# as they should all be single peptides since they are not grouped around mRNA translating
# Returns a dictionary with image names as keys and values as median of single peptide in that image
def protAvg(mrna_spots, protein_spots, threshold_distance, num_prot_avg, min_prot_avg):
    # Initialize empty list which will be filled with all image names and used to check if each image has a median protein value
    list_of_all_image_names = []
    # Initialize empty dictionary to be filled with a list of protein amplitudes not within an mRNA by image name (key)
    dict_of_prot_brightness = {}
    # Initialize empty dictionary to be filled with median single protein amplitudes by image name (key)
    dict_of_brightness_median = {}
    
    for protein_entry in protein_spots:
        # If image name that protein belongs to is in dictionary and has enough points to take median for skip over protein
        if protein_entry[0] in dict_of_prot_brightness:
            if dict_of_prot_brightness[protein_entry[0]][0] > num_prot_avg:
                continue
        # If image name has not already been added to list_of_all_image_names add it to the list
        if not protein_entry[0] in list_of_all_image_names:
            list_of_all_image_names += [protein_entry[0]]
        # Extract protein coordinates from the data. They are 3rd, 4th and 5th rows
        protein_coordinates = (protein_entry[2:5].astype(float))
        within_mrna = False
        for mrna_entry in mrna_spots:
            # Check and make sure that mRNA and protein are in the same image and the same cell number
            if mrna_entry[0].split(mRNA_channel_path)[0] == protein_entry[0].split(protein_channel_path)[0] and mrna_entry[1] == protein_entry[1]:
                mrna_coordinates = (mrna_entry[2:5].astype(float))
                # Calculate distance between mRNA and protein
                distance = (math.sqrt(np.sum((mrna_coordinates - protein_coordinates) ** 2)))
                # If distance is less than the threshhold for translation update the within mRNA variable to True and exit loop
                if distance <= int(threshold_distance):
                    within_mrna = True
                    break  
        # If protein is not translating add it to list for specified image name within the dictionary dict_of_prot_brightness
        # Note: dictionary = [counter, [list of prot]]. The counter here is to explicitly count the number of proteins in list
        if not within_mrna:
            if protein_entry[0] in dict_of_prot_brightness:
                dict_of_prot_brightness[protein_entry[0]][1] += [float(protein_entry[5])]
                dict_of_prot_brightness[protein_entry[0]][0] += 1
            # If this is the first protein for the image create a new key and add protein
            else:
                dict_of_prot_brightness[protein_entry[0]] = [1,[float(protein_entry[5])]]
    # Loop through the image names used as keys in the dictionary and take the median of the protein amplitudes if there is enough 
    # proteins to take a median (number of proteins >= min_prot_avg)
    for image in dict_of_prot_brightness:
        if dict_of_prot_brightness[image][0] >= min_prot_avg:
            dict_of_brightness_median[image] = np.median(dict_of_prot_brightness[image][1])
    # Add image names that did not have enough proteins to get median to dictionary with a value of -1, signaling that this file must be skipped
    for image_name in list_of_all_image_names:
        if not image_name in dict_of_brightness_median:
            dict_of_brightness_median[image_name] = -1
    return dict_of_brightness_median

def treatment_data(treatment_list, data, outputFileName):
    # If the user elects to use the default treatment settings by entering zero set treatment_list to default
    if treatment_list[0] == "0":
        treatment_list[0] == ["Ctrl", "MG132"]
    
    # Initialize two empty numpy arrays to be propegated by our data
    new_count_data = np.full((len(data), len(treatment_list)*21),None, dtype = object)
    new_translate_data = np.full((len(data),len(treatment_list)*2),None, dtype=object)
    # Initialize two arrays of zeros of length of the number of different treatments there are
    counter_list = [0]*len(treatment_list)
    counter_skip = [0]*len(treatment_list)
    # Generate file names of the csv files we will write
    count_filepath_csv = str(os.getcwd()+"/"+outputFileName+"_control_heat_count_by_image.csv")
    translating_filepath_csv = str(os.getcwd()+"/"+outputFileName+"_control_heat_translating_by_image.csv")
    
    # Initialize empty header lists
    headers_count = []
    headers_translate = []
    # Use loop to create the list of headers for count and translating csv, final format of headers can be seen in sample files
    for i in range(1,22):
        # Loop through treatment list
        for j in range(len(treatment_list)):
            if i == 21:
                headers_count += ["20+ " + treatment_list[j]]
            else:
                headers_count += [str(i) + " " + treatment_list[j]]
            if i == 1:
                headers_translate += ["Translating " + treatment_list[j]]
            if i == 2:
                headers_translate += ["Not Translating " + treatment_list[j]]

    # Loop through treatments and the data to sort data by treatment 
    for i in range(len(treatment_list)):
        for row in data:
            # If treatment key word is in the name of the file
            if treatment_list[i] in row[0]:
                # Get total mRNA to use to determine percentage of untranslating and translating
                total_mrna = float(row[1])+float(row[2])
                # Add translating percentage to the numpy array of translation data, row determined by the amount of that treatment has
                # already been done Coloumn is determined by location in treatment list, see example spreadsheet for more clarity
                new_translate_data[counter_list[i], i] = float(row[2])/total_mrna*100
                # Adds non-translating percentage to the translation data numpy array
                new_translate_data[counter_list[i], len(treatment_list) + i] = float(row[1])/total_mrna*100
                # Ignores row of data if there are no translating mRNA, (Can't get percentages by number of translating mRNA)
                if int(row[2]) != 0:
                    # Loop through coloumns in row of data
                    for j in range(0,21):
                        # Add mRNA count data to new numpy array
                        new_count_data[counter_list[i]-counter_skip[i], len(treatment_list)*j+i] = float(row[j+4])/float(row[2])*100
                # If row is skipped over add one to skip counter of that treatment, this ensures that data printed with no empty rows
                else:
                    counter_skip[i] += 1
                # Increment counter to go to next row when adding data for this specfic treatment
                counter_list[i] += 1
    # Write newly orginized data into csv files, with the headers we generated and the name decided
    write_csv(translating_filepath_csv, new_translate_data, headers_translate)
    write_csv(count_filepath_csv, new_count_data, headers_count)

# Take user input, a complete guide to user input is explained in the instructions document for this code
channel1Path = input("What is the full path of your mRNA channel? ") #mRNA channel output file from fishquant - summary file all spots 
channel2Path = input("What is the full path of your peptide channel? ") #peptide channel output file from fishquant - before threshold cutoff
# Protein and mRNA channel paths are used to remove them from respective image names so it can be checked that protein and mRNA image
# names are the same when comparing coordinates
protein_channel_path = input("What is the protein channel path (ex. '640S.tif'): ")
mRNA_channel_path = input("What is the mRNA channel path (ex. '555S.tif'): ")
# A list of treatment names is used to determine what image belongs to which treatment so they can be orginized nicely in spreadsheets
treatments = input("Input treatment labels seperated by commas without spaces. To use default input 0: ").split(",")


try:
    if os.path.exists(channel1Path):

        os.chmod(channel1Path, 0o666)
        print("File permissions modified successfully!")
    else:
        print("File not found:", channel1Path)
except PermissionError:
    print("Permission denied: You don't have the necessary permissions to change the permissions of this file.")
    
    
# Load mRNA and protein channels in as ch1 and ch2 respectively
ch1=np.loadtxt(channel1Path, skiprows=14, dtype=str)   
ch2=np.loadtxt(channel2Path, skiprows=14, dtype=str)   

knowMean = input("What is the mean of peptide? If you don't know enter 0, if you want to use different means for different images enter -1, otherwise enter the mean for all images: ")
# User has elected to input a list of image names and an associated list of mean protein amplitudes
if knowMean == "-1":
    mean = {}
    while True:
        # Ask user for list of image names and list of mean values, if they are not same length as them to try again
        try:
            list_of_image_names_initial = input("List of image names seperatd by commas: ")
            # Rid user input of single quotes and spaces, this ensures compatibility with printed output, allowing user to copy and paste
            # split input by commas into a list
            list_of_image_names = list_of_image_names_initial.replace("'", "").replace(" ", "").split(",")
            # Similiar to image names but only spaces have to be removed
            list_of_mean_values = input("List of mean values for file names in the same orderas the image names above: ").replace(" ", "").split(",")
            # Raise exception if list are not same length
            if len(list_of_image_names) != len(list_of_mean_values):
                raise Exception("The length of your image name list and mean list are of different lengths, you must have a mean value for each image name")
            # break out of loop if exception is not reached
            break
        except Exception as e:
            continue
    # Propogate dictionary with image names as keys and coresponding protein amplitude as the value
    for i in range(len(list_of_image_names)):
        mean[list_of_image_names[i]] = float(list_of_mean_values[i])

thresholdedDistance = input("what is the distance threshold between mRNA and peptide?")   #if the distance b/w peptide and mRNA is less than or equal to 100
outputFileName = input("Save the output file as:")

# If user does not input symbolic value but rather an actual positive value mean, propogate the dictionary with image names as keys and
# the protein amplitude as the mean for all images, less specific and likely less acurate
if (knowMean!="0" and knowMean != "-1"):
    mean = {}
    for protein_entry in ch2:
        if protein_entry[0] not in mean: 
            mean[protein_entry[0]] = float(knowMean)

# User inputs zero indicating they do not know what the mean is and need to calculate it from the data
elif (knowMean == "0"):
    definition_of_single = input("How would you like to define a single, gaussian (type g) or mrna away from protien (type p): ")
    if (definition_of_single == "g"):
        thresholdedSpotFile = input("What is the full path of the summary thresholded spots? ")
        mean_estimation_file=np.loadtxt(thresholdedSpotFile, skiprows=14,dtype=str) #peptide channel output file from fishquant for SINGLES
        mean=np.mean(np.array(mean_estimation_file[:,5].astype(float)))
    elif (definition_of_single == "p"):
        num_prot_avg = int(input("How many proteins do you want to take the median of to get the amplitude of a single protein: "))
        min_prot_avg = int(input("What is the minimum number of points that need to be not translating for image to be included in analysis (minimum number of points averaged to get single peptide amplitude): "))
        mean = protAvg(ch1,ch2,thresholdedDistance, num_prot_avg, min_prot_avg)
    else:
        raise ValueError("You must input either 'g' for gaussian or 'p' for mrna away from protien")

# Print out image names and mean for user reference, they can use these lists as input in future program runs to save
# time spen on median calculations
key_list = []
mean_list = []
for key in mean:
    print("Mean peptide intensity of image", key, "is:", mean[key])
    key_list += [key]
    mean_list += [mean[key]]
print("Image name list and mean value list:")
print("List of image names: ", key_list)
print("List of mean values:", mean_list)

# Open txt file that results will be written to
file=open(outputFileName+"_fishQuantTransltion_Results.txt",'w+')
# Write headers to txt file
file.write('\t'.join(['image_name','cell_number','mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance','X Coord', "Y Coord", "Z Coord"])+"\n")
# Create header list for csv version of this file
headers = ['image_name','cell_number','mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance', 'X Coord', "Y Coord", "Z Coord"]
# Initialize list which will contain lists repersenting rows of the csv file
summary_total = []
# Loop through the length of the mRNA and protein entries
for entry in range(ch1.shape[0]):
    for peptide_entry in range(ch2.shape[0]):
        # File name is the first coloumn of a given row
        file_name = ch2[peptide_entry,0]
        # Only use data if a acceptable peptide amplitude median exists for this given file name
        if mean[file_name] > 0:
            # Check to ensure that the protein and mRNA molecules are within the same cell and within the same image
            # Split function removes the channel pathway so image names can be compared
            if (ch1[entry,1]==ch2[peptide_entry,1]) and (ch1[entry,0].split(mRNA_channel_path)[0]==ch2[peptide_entry,0].split(protein_channel_path)[0]):
                # Get the x y z coordinates from the files (located in 2,3 and 4th column)
                p1= ch1[entry,2:5].astype(float)
                p2= ch2[peptide_entry,2:5].astype(float)
                # Extract image name and cell number
                image_name=ch1[entry,0].astype(str)
                cell = ch1[entry,1].astype(str)
                # Calculate distance between mRNA and protein
                distance=math.sqrt(((p1[0]-p2[0]))**2+((p1[1]-p2[1]))**2+((p1[2]-p1[2]))**2)
                # If the distance is underneath the user generated threshold for translating, add new row with data
                if distance<=float(thresholdedDistance):
                    summary=[]
                    summary.append(str(image_name))
                    summary.append(str(cell))
                    summary.append(str(entry))
                    summary.append(str(peptide_entry))
                    # Peptide amplitude
                    summary.append(str(ch2[peptide_entry,5]))
                    # Normalized peptide amplitude
                    summary.append(str(ch2[peptide_entry,5].astype(float)/mean[file_name]))
                    summary.append(str(distance))
                    # X coordinate
                    summary.append(str(ch2[peptide_entry,3]))
                    # Y coordinate
                    summary.append(str(ch2[peptide_entry,2]))
                    # Z coordinate
                    summary.append(str(ch2[peptide_entry,4]))
                    # Add summary list to total summary list to be turned into a csv
                    summary_total += [summary]
                    # Write row to txt file
                    file.write("\t".join(summary)+"\n")
file.close()
# Create csv file name and call write_csv() to make new csv file with data
outputfilepath_csv = str(os.getcwd()+"/"+outputFileName+"_fishQuantTransltion_Results.csv") 
write_csv(outputfilepath_csv,summary_total,headers)


################################remove repeated mRNA, based the brightness ################################
# Open file we just created to narrow down the protein dots associated with mRNA to one, only keeping the one that is the brightest (only one protein dot per mRNA)
outputfilepath = str(os.getcwd()+"/"+outputFileName+"_fishQuantTransltion_Results.txt")
summary=np.loadtxt(outputfilepath, skiprows=1, dtype=str)
# Create a list of unique mRNA entry numbers
uniquemRNAs=np.unique(summary[:,2])
# Open new txt file and write to it headers
textAmp=open("tmp_amplitude.txt","w+")
textAmp.write("\t".join(['image_name','cell_number','mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance'])+"\n")

# Loop through unique mRNA entry numbers
for mRNA in uniquemRNAs:
    # Create a list of roww indexs that have this specific mRNA entry number
    repeatedEntries=np.where((summary[:,2]==mRNA))
    # If more than one proetin associated with a single mRNA only keep the one that has the highest amplitude   
    if len(repeatedEntries[0])>1:
        # Initialize entry_id (row number) and max amp value to -1
        entry_id = -1
        max_pep_amp = -1
        # Loop through indexs that have this specific mRNA entry number
        for entry in repeatedEntries[0]:
            # If this row has a higher peptide amplitude keep index of row and update max_pep_amp to the row's amplitude
            if float(summary[entry,5]) > max_pep_amp:
                entry_id = entry
                max_pep_amp = float(summary[entry,5])
        # Write the row for this specific mRNA entry number that has the highest peptide amplitude
        textAmp.write("\t".join(summary[entry_id,:]))
        textAmp.write("\n")
    # If there is only 1 row for the mRNA entry number write that row to the txt file
    else:
        entry_id=repeatedEntries[0][0]
        textAmp.write("\t".join(summary[entry_id,:]))
        textAmp.write("\n")
# Close them txt file
textAmp.close()


###############################remove repeated peptides, based on distance ###############################
# In this section we will get rid of rows with repeated peptide entries by deleting all but the one where the peptide is closest to it's mRNA
# Load in temp file created above
summaryAmpDist=np.loadtxt("tmp_amplitude.txt", skiprows=1, dtype=str)
# Create list of unique peptide entries
uniquemPEPs=np.unique(summaryAmpDist[:,3])
# Open new txt file to write results to and write headers to the file
textDistance=open("tmp_distance.txt","w+")
textDistance.write("\t".join(['image_name','cell_number', 'mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance'])+"\n")

# Loop through the list of unique peptide entries
for PEP in uniquemPEPs:
    # Create a list of indexs where peptide entry is PEP
    repeatedEntries=np.where((summaryAmpDist[:,3]==PEP))    
    # If multiple peptide have this specfic peptide entry number determine which mRNA it is closest to
    if len(repeatedEntries[0])>1:
        # Initial entry id to -1 and minimum distance to the max value of an integer
        entry_id = -1
        distance_min = sys.maxsize
        # Loop through the repeated indexs of specific peptide
        for entry in repeatedEntries[0]:
            # If distance of row is lower then distance_min update entry_id to this row and update distance_min
            if float(summaryAmpDist[entry, 6]) < distance_min:
                entry_id = entry
                distance_min = float(summaryAmpDist[entry,6])
        # Write row with lowest distance to the txt file
        textDistance.write("\t".join(summaryAmpDist[entry_id,:]))
        textDistance.write("\n")
    # If there is only a single row with the specific peptide entry number write row to txt file
    else:
        entry_id=repeatedEntries[0][0]
        textDistance.write("\t".join(summaryAmpDist[entry_id,:]))
        textDistance.write("\n")
textDistance.close()


############################### for prism output --> translation efficiency ###############################
# This section creates a txt and csv file containing the number of translating and non translating mRNA by cell
# Load in txt file created in the previous section, these rows will only contain spots that are within the translating threshold
summaryAmpDist=np.loadtxt("tmp_distance.txt", skiprows=1, dtype=str)
# Concatenate the image name and the cell number together to get a unique image/cell name, from txt file created previously
unified_translating_cellnames=np.char.add(summaryAmpDist[:,0].astype(np.str_),summaryAmpDist[:,1].astype(np.str_))
# Do same concatenation as above but from orginal txt file containing mRNA information
unified_total_cellnames=np.char.add(ch1[:,0].astype(np.str_),ch1[:,1].astype(np.str_))
# Create txt and csv file path names
translatingFilepath = str(os.getcwd()+"/"+outputFileName+"_cell_translatingmRNA.txt")
translatingFilepath_csv = str(os.getcwd()+"/"+outputFileName+"_cell_translatingmRNA.csv")
# Open txt file and write header names
text=open(translatingFilepath,"w+")
text.write("\t".join(['image_name','cell_number','translating','untranslating'])+"\n")
# Create header list for csv file output
headers = ['image_name','cell_number','translating','untranslating']
# Initialize empty data list to hold rows of csv file as lists
translation_data = []

# Loop through the unique file_names in unified file name list
for name in np.unique(unified_translating_cellnames):
    # The length of the list where the cell name and image name are the same repersents the number of translating proteins in that cell
    translating_count=len(np.where(unified_translating_cellnames==name)[0])
    # The length of the mRNA channel will have the toal number of mRNA with the same cell name and image name
    total_count=len(np.where(unified_total_cellnames==name)[0])
    # Total number of mRNA with this cell/image name minus the number of translating mRNA will give the untranslating mRNA
    untranslating_count = total_count - translating_count
    # Split word at .tif back into image name and cell number
    itms=name.split('.tif')
    # Write to csv
    text.write('\t'.join([itms[0]+'.tif',itms[1]])+"\t"+str(translating_count)+"\t"+str(untranslating_count)+"\n")
    file_name = itms[0]+'.tif'
    # Add row to data to be eventually written into csv
    translation_data += [[file_name, itms[1], translating_count, untranslating_count]]
# Close txt file and write csv with file name, data, and headers created above
text.close()
write_csv(translatingFilepath_csv, translation_data, headers)

############################### for prism output --> ribosome occupency ###############################
# This section generates the various different count spreadsheets as well as the one count txt file
# Open txt file created in the previous section to extract data from
summaryAmpDist=np.loadtxt("tmp_distance.txt", skiprows=1, dtype=str)
# Get unified cell names so that rows in the same image and cell can be identified
unified_cellnames=np.char.add(summaryAmpDist[:,0].astype(np.str_),summaryAmpDist[:,1].astype(np.str_))
# Create a list of image names, to be used during csv creation
all_file_names = summaryAmpDist[:,0]
# Initialize an empty array where csv file data will be held for translating and untransalting by image name rather then cell and image name
data_combined_by_file_name = []

# The purpose of this loop is to combine translating and untranslating mRNA counts for all cells in an image into total translating and
# untranslating for the whole image
for row in translation_data:
    # If data_combined_by_file_name array is empty add a row with image name, translating mRNA, and untranslating mRNA
    if not data_combined_by_file_name:
        data_combined_by_file_name += [[row[0], int(row[2]), int(row[3])]]
    # If array has been initialized
    else:
        # If image name has not already been written into array add row without cell number from translating data array
        if not row[0] in np.array(data_combined_by_file_name)[:,0]:
            data_combined_by_file_name += [[row[0], int(row[2]), int(row[3])]]
        # If there already exists a row in data_combined_by_file_name with this rows image name, add translating and untranslating data
        else:
            # Loop through the rows in new array
            for row_index in range(len(data_combined_by_file_name)):
                # Find row where the image name is the same
                if data_combined_by_file_name[row_index][0] == row[0]:
                    # Add translating and untranslating mRNA from row in translating_data to row with same image name in new spreadsheet
                    data_combined_by_file_name[row_index][1] += int(row[2])
                    data_combined_by_file_name[row_index][2] += int(row[3])

# Create file paths for txt file and by image and by cell count csv files and by image and by cell percentage count csv files
countfilepath = str(os.getcwd()+"/"+outputFileName+"_count.txt")
countfilepath_by_cell_csv = str(os.getcwd()+"/"+outputFileName+"_count_by_cell.csv")
countfilepath_by_image_csv = str(os.getcwd()+"/"+outputFileName+"_count_by_image.csv")
countfilepath_by_cell_percent_csv = str(os.getcwd()+"/"+outputFileName+"_count_by_cell_percent.csv")
countfilepath_by_image_percent_csv = str(os.getcwd()+"/"+outputFileName+"_count_by_image_percent.csv")
# Open txt file and write the headers to teh file
text=open(countfilepath,"w+")
text.write('\t'.join(['image_name','cell_number','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','20+'])+"\n")
# Create the header list for both csv files
headers_by_cell = ['image_name','cell_number','Not translating', 'Total translating','Median single protien amp', '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','20+']
headers_by_image = ['image_name','Not translating', 'Total translating','Median single protien amp', '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','20+']
# Initialize empty arrays for all our csv data to be created
new_data_by_cell = []
new_data_by_image =[]
percentage_count_data_by_cell = []
percentage_count_data_by_image = []
# Create an empty list for file names that have been used
file_name_used = []
# Initialize counters for image and cell csv files
counter_for_by_image = 0
counter_for_by_cell = 0
# Loop through all unique unified cell/image names
for name in np.unique(unified_cellnames):
    # Split at the .tif mark to break unified cell name into image and cell names
    itms=name.split('.tif')
    # Add .tif back to get complete image name
    file_name = itms[0]+'.tif'
    # Change the channel pathway on the image name to the protein channel pathway as that is what acts as a key for median single amplitude peptide values
    key = name.split(mRNA_channel_path)[0]
    key += protein_channel_path
    # If the image does not have a valid single peptide amplitude (-1) ignore this image
    if mean[key] > 0:
        # Initialize empty lists that will be used to create rows of the csv files
        by_cell_counts =[]
        by_image_counts = []
        # If the cell name has been used skip over
        if not file_name in file_name_used:
            # Check for which indexes in previously created translating mRNA file, have image name equal to the image name that we are
            # currently looping through and have a normalized peptide amplitude less than 0.5, the length of this index list will be the
            # number of translating proteins in the image with peptide ampltude of under 0.5
            by_image_counts += [len(np.where((all_file_names==file_name) & (summaryAmpDist[:,5].astype(float)<0.5))[0])]
            # Similiar to above with different normalized peptide amplitude range values
            by_image_counts += [len(np.where((all_file_names==file_name) & (summaryAmpDist[:,5].astype(float)>=0.5) & (summaryAmpDist[:,5].astype(float)<=1.0))[0])]
            by_image_counts += [len(np.where((all_file_names==file_name) & (summaryAmpDist[:,5].astype(float)>1.0) & (summaryAmpDist[:,5].astype(float)<=(0.74*2)))[0])]
            # Loop through to complete the rest of the different range values
            for i in range (2,21):
                by_image_counts += [len(np.where((all_file_names==file_name) & (summaryAmpDist[:,5].astype(float)>(0.74*i)) & (summaryAmpDist[:,5].astype(float)<=(0.74*(i+1))))[0])]
            # Add last set of values for anything over a set point
            by_image_counts += [len(np.where((all_file_names==file_name) & (summaryAmpDist[:,5].astype(float)>(0.74*20)))[0])]
            # Add file name that was just completed to list of completed file names as to not repeat any
            file_name_used += [file_name]
            # Add row to data by image array, first coloumn file name, second coloumn non-translating mRNA (taken from previously created
            # spreadsheet but with the addition of weak amplitude peptides), third coloumn translating mRNA (taken from previously created
            # spreadsheet but subtracting the weak amplitude peptides), remaining coloumns are number of peptide spots at each size (1-20+)
            new_data_by_image += [[file_name, data_combined_by_file_name[counter_for_by_image][2]+by_image_counts[0], data_combined_by_file_name[counter_for_by_image][1]-by_image_counts[0], mean[key], by_image_counts[1], by_image_counts[2], by_image_counts[3], by_image_counts[4], by_image_counts[5], by_image_counts[6], by_image_counts[7], by_image_counts[8], by_image_counts[9], by_image_counts[10], by_image_counts[11], by_image_counts[12], by_image_counts[13], by_image_counts[14], by_image_counts[15], by_image_counts[16], by_image_counts[17], by_image_counts[18], by_image_counts[19], by_image_counts[20], by_image_counts[21]]]
            # Create similiar row using one we already made, but with percentages instead of counts of peptide spots at each size
            percentages_by_image = []
            # Check to ensure number of translating mRNA is above zero before trying to divide by it, divide by zero error
            if new_data_by_image[counter_for_by_image][2] > 0:
                for value in new_data_by_image[counter_for_by_image][4:]:
                    percentages_by_image += [value/new_data_by_image[counter_for_by_image][2]*100]
                percentage_count_data_by_image += [new_data_by_image[counter_for_by_image][0:4] + percentages_by_image]
            # If number of translating mRNA is 0 then row in perecent spreadhseet is the same as row in count spreadsheet
            else:
                percentage_count_data_by_image = [new_data_by_image[counter_for_by_image]]
            # Increment count for by image csv array to make sure the correct row is written in
            counter_for_by_image += 1

        # Very similiar process to what was done by image above but this time unified cellname is used to do it by image
        by_cell_counts += [len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)<0.5))[0])]
        by_cell_counts += [len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>=0.5) & (summaryAmpDist[:,5].astype(float)<=1.0))[0])]
        by_cell_counts += [len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>1.0) & (summaryAmpDist[:,5].astype(float)<=(0.74*2)))[0])]
        for i in range (2,20):
            by_cell_counts += [len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>(0.74*i)) & (summaryAmpDist[:,5].astype(float)<=(0.74*(i+1))))[0])]
        by_cell_counts += [len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>(0.74*20)))[0])]
        # One difference to above is the second coloumn is teh cell number
        new_data_by_cell += [[file_name, itms[1], translation_data[counter_for_by_cell][3]+by_cell_counts[0], translation_data[counter_for_by_cell][2]-by_cell_counts[0], mean[key], by_cell_counts[1], by_cell_counts[2], by_cell_counts[3], by_cell_counts[4], by_cell_counts[5], by_cell_counts[6], by_cell_counts[7], by_cell_counts[8], by_cell_counts[9], by_cell_counts[10], by_cell_counts[11], by_cell_counts[12], by_cell_counts[13], by_cell_counts[14], by_cell_counts[15], by_cell_counts[16], by_cell_counts[17], by_cell_counts[18], by_cell_counts[19], by_cell_counts[20], by_cell_counts[21]]]
        percentages_by_cell = []
        if new_data_by_cell[counter_for_by_cell][3] > 0:
            for value in new_data_by_cell[counter_for_by_cell][5:]:
                percentages_by_cell += [value/new_data_by_cell[counter_for_by_cell][3]*100]
            percentage_count_data_by_cell += [new_data_by_cell[counter_for_by_cell][0:5] + percentages_by_cell]
        else:
            percentage_count_data_by_cell += [new_data_by_cell[counter_for_by_cell]]
        counter_for_by_cell += 1
        # Write by cell results into txt file
        text.write("\t".join([itms[0]+'.tif', itms[1]]) + "\t" + \
           "\t".join(str(by_cell_counts[i]) for i in range(1, len(by_cell_counts))) + \
           "\n")
        
# Write all results created above (by_cell, by_image, by_cell_percent, by_image_percent) into csv files
write_csv(countfilepath_by_cell_csv, new_data_by_cell, headers_by_cell)
write_csv(countfilepath_by_image_csv, new_data_by_image, headers_by_image)
write_csv(countfilepath_by_cell_percent_csv, percentage_count_data_by_cell, headers_by_cell)
write_csv(countfilepath_by_image_percent_csv, percentage_count_data_by_image, headers_by_image)

# Close txt file
text.close()

# Call treatment_data function described above to create new spreadsheet that is better orginized for analysis 
treatment_data(treatments, new_data_by_image, outputFileName)

############################### clean up ###############################
os.remove('tmp_amplitude.txt')
os.remove('tmp_distance.txt')

print("done")