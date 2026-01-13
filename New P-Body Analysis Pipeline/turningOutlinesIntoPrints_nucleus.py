# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:49:31 2024

@author: yorub
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_fill_holes
from scipy.interpolate import interp1d
import os
from mahotas.polygon import fill_polygon

## folder containing outline files

#folder = r"C:\Users\yorub\Downloads\outline"
#skel_folder_name = r"C:\Users\yorub\Downloads\SkeletonsAndPrints"
skel_folder_name = input("Root folder \n")
skel_folder_name = os.path.join(skel_folder_name, "SkeletonAndPrint")
folder = input("Input folder with outline files \n")

files = os.listdir(folder)

if not os.path.exists(skel_folder_name):
    os.makedirs(skel_folder_name)

for file in files:
    path = os.path.join(folder, file)
    filename = file[:file.index("__outline")]
    
    with open(path) as f:
        lines = f.read()
        cells = lines.split("Nucleus_START")[1:]
        
        cell_num = 1
        for cell in cells:
            lines_in_cell = cell[:cell.index("Nucleus_END")]
            lines_yPOS = lines_in_cell[lines_in_cell.index("X_POS"):lines_in_cell.index("Y_POS")]
            lines_xPOS = lines_in_cell[lines_in_cell.index("Y_POS"):lines_in_cell.index("Z_POS")]
            
            lines_xPOS = lines_xPOS.split()[1:]
            lines_xPOS += [lines_xPOS[0]]
            lines_yPOS = lines_yPOS.split()[1:]
            lines_yPOS += [lines_yPOS[0]]
            
            cell_xpositions = np.array(lines_xPOS, dtype=int)
            cell_ypositions = np.array(lines_yPOS, dtype=int)
            
            mask = np.zeros((2304,2304))
            
            fill_polygon( np.array([cell_xpositions, cell_ypositions]).T, mask )
            '''
            mask[cell_xpositions, cell_ypositions] = 1
            for i in range(len(cell_xpositions)-1):
                xcurr, ycurr = cell_xpositions[i], cell_ypositions[i]
                xnext, ynext = cell_xpositions[i+1], cell_ypositions[i+1]
                
                func = interp1d([xcurr, xnext], [ycurr, ynext])
                all_xs = np.arange(xcurr, xnext)
                all_ys = func(all_xs)
                mask[all_xs, np.array(all_ys, dtype=int)] = 1
                
            mask = binary_fill_holes(mask)
            '''
            plt.imsave(os.path.join(skel_folder_name, filename + f"_nuc_{cell_num}.gif"), mask)
            cell_num += 1
            