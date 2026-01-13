# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 12:17:07 2022

@author: Administrator
"""
import os
import re
#import bigfish.multistack as multistack
import pandas as pd
#import bigfish.stack as stack




def getImagesDir():
    # Getting the image folder
    print('Input path to working directory.')
    print("\n Full path to working directory: ", end = "")
    directory = input()
      
    while not os.path.isdir(directory):
      print("\n Incorrect path, please try again:")
      print("\n Full path to folder: ", end = "")
      directory = input()
    return directory

def regexMatchFilter(lis, pattern, show=False):
    '''
    Parameters
    ----------
    lis : list
        list of filenames.
    pattern : Regex match object
        Naming format to match filenames to

    Returns
    -------
    matches : list
        List of filenames that match the format.
    '''
    matches = []
    i=0
    for l in lis:
      
      m = re.fullmatch(pattern, l)
      if m: 
          matches.append(m)
          if show:
              print("Index:", i, ". Filename: ",m.group(0)[:-4])
          i+=1
    return matches

def getMaxProjectionImages(directory, return_chans = False):
    '''
    Function to get match object based on the file naming format to segment multiple images at once
    Asks user for image naming structure and returns 2D list of all images: first list is of all DAPI images and
    second is of all cytoplasm stained images
    Parameters
    ----------
    directory : str
        Full path to folder containing MAX projection images.

    Returns
    -------
    list
        Two dimensional list: first list is of all DAPI images and
        second is of all cytoplasm stained images
    '''
    #for now we just input one image each
    lis = os.listdir(directory)
    
    dapi = input("Enter name of DAPI channel (stained nucleus):\n")
    #dapi_pattern = re.compile(".*MAX_(.*?)_(.*?)_(.*?)_.*_(.*?)_(.*?)_(?:xy|XY)([0-9]?[0-9])_" + dapi + ".gif")
    dapi_pattern = re.compile(".*MAX_.*" + dapi + ".gif")
    dapi_files = regexMatchFilter(lis, dapi_pattern)
    
    
    map2 = input("Enter name of RNA channel: \n")
    #map2_pattern = re.compile(".*MAX_(.*?)_(.*?)_(.*?)_.*_(.*?)_(.*?)_(?:xy|XY)([0-9]?[0-9])_" + map2 + ".gif")
    map2_pattern = re.compile(".*MAX_.*" + map2 + ".gif")
    map2_files = regexMatchFilter(lis, map2_pattern)
    
    if return_chans:
        return [dapi_files, map2_files, dapi, map2]
    return [dapi_files,map2_files] #DAPI, cytoplasm


def getSegmentationMaxProjections(directory):
    lis = os.listdir(directory)
    
    
    dapi = input("Enter name of DAPI channel")
    #dapi = "395S"
    dapi_pattern = re.compile(".*MAX_.*" + dapi + ".gif")
    dapi_files = regexMatchFilter(lis, dapi_pattern)
    
    channels = []
    chan = input("Enter name of a channel which shows the cytoplasm. If no more, press 0.")
    
    while chan != "0":
        channels.append(chan)
        chan = input("Enter name of a channel which shows the cytoplasm. If no more, press 0.")
    
    #channels = ["740S","640S","555S"]
    
    cytoplasm_images = []
    for chan in channels:
        print(len(cytoplasm_images))
        chan_pattern = re.compile(".*MAX_.*" + chan + ".gif")
        chan_files = regexMatchFilter(lis, chan_pattern)
        print("chanfiles", len(chan_files))
        cytoplasm_images.append(chan_files)
    
    return [dapi_files, cytoplasm_images]

def get_TIF_and_MAX_Images(directory, message, show=False, return_chan=False):
    '''gets max projection and TIF files
    message should ask user to input correct channel name
    str, str -> 2D list of match objects'''
    
    lis = os.listdir(directory)
    
    chan = input(message)
    
    max_pattern = re.compile(".*MAX_.*" + chan + ".gif")
    max_files = regexMatchFilter(lis, max_pattern, show=False)
    
    
    tif_pattern = re.compile(".*_" + chan + ".tif")
    tif_files = regexMatchFilter(lis, tif_pattern, show=show)
    
    if return_chan:
        return max_files,tif_files, chan
    else:
        return [max_files,tif_files] 

def getMaxImagesAndFilenames(chan1, chan2):
    '''the var naming assumes chan1 is the MAP2 stain and chan2 is rna stain, but this would generall work for any 2 channels
    Returns [maxprojection images of chan1, TIF files of chan1, chan1 filenames, maxprojection images of chan2, TIF files of chan2, chan2 filenames]'''
    map2, rna = chan1, chan2 #a list of Match objects for every map2 and RNA max projection image
    dend_mips = [] #list of ndarrays, max projections of dendrite/MAP2 stains
    dends = [] #list of ndarrays, TIF files of dendrite/MAP2 stains
    dend_filenames, rna_filenames = [], [] #list of names of files excluding file extension and excluding MAX_ (ie only contains all the experiment info)


    for match in map2:
        dend_mips.append(stack.read_image(os.path.join(path_input, match.group(0))))
        
        filename = match.group(0)[4:-4] #Exclude the beginning "MAX_" and the ending ".gif"
        dend_filenames.append(filename)
        
        dends.append(stack.read_image(os.path.join(path_input,filename+".tif")))
        
    rna_mips = [] #list of max projections of rna stain
    rnas = [] #list of TIF files of rna stains
    
    
    for match in rna:

        rna_mips.append(stack.read_image(os.path.join(path_input, match.group(0))))
        
        filename = match.group(0)[4:-4]
        rna_filenames.append(filename)

        rnas.append(stack.read_image(os.path.join(path_input,filename+".tif")))
        
    if len(rna_mips) != len(dend_mips):
        print("There's something wrong here!")
        
    return [dend_mips, dends, dend_filenames, rna_mips, rnas, rna_filenames]

def getYesorNo(message):
    print(message)
    ans = input()
    possible_ans = ['Y','y','n','N', '0','1']
    while ans not in possible_ans:
        print("You did not input Y or N. Please try again.")
        ans = input(message)
    
    pos_ans = ['Y','y','1']
    neg_ans = ['N','n','0']
    
    if ans in pos_ans:
        return True
    if ans in neg_ans:
        return False

def saveDfDictToExcel(data_dict, name):
  with pd.ExcelWriter(name) as writer:
    for div_treat, df in data_dict.items():
      cols = df.columns.values.tolist()
      df.to_excel(writer, sheet_name=div_treat, columns=cols, index=False)