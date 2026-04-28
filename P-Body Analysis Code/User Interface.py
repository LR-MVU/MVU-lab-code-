from __future__ import annotations

import os
from pathlib import Path

from pbody_charecteristics_2D import PBody_Charecteristics
from pbody_mrna_metrics_2D import run_pbody_mrna_metrics

#Need to alter cropping to ensure that cell numbers are of fixed length, 010 vs 0010

# Edit these values directly before running.
#path_root = Path(r"Z:\Organized Codes\Re-organized P-Body Code\Test Folder\bigFISH_In")
#conditions = ["Dox-_2HR", "Dox-_CT", "Dox+_CT"]
path_root = Path(r"Z:\Organized Codes\Re-organized P-Body Code\Test Folder\bigFISH_In")
#conditions = ['Puro_HS','Puro_2HR', 'Puro_3HR', 'Puro_4HR', 'Puro_6HR', 'Puro_8HR', 'Puro_10HR', 'Puro_CT', 'WT_2HR', 'WT_3HR', 'WT_4HR', 'WT_6HR', 'WT_8HR', 'WT_10HR', 'WT_CT', 'WT_HS']
conditions = ["Dox-_2HR", "Dox-_CT", "Dox+_CT"]
#['Ctrl', 'HS', '2H_R', '3H_R', '4H_R', '6H_R', '8H_R', '10H_R', '12H_R', '24H_R']

diffraction = 0  #140
pixel_res = 107.5
pbody_chan = '640S'

rna_chan = '555S'
spot_source = 'declustered' # Takes values "declustered" or "raw"
#Formatting for RNA csv files must match the following
# Declustered: f"{mask_basename}_mRNA{rna_chan}_spots.csv"
# Raw: f"{mask_basename}_mRNA{rna_chan}_spots_ind.csv"
# If your excel file naming does not follow this format, you can change it manually ,
# On line 108 or 110 of the "pbody_mrna_metrics_2D.py" code file, make sure to revert afterwards

cell_mask_chan = 'print'
nuc_mask_chan = 'nuc'
pb_folder = 'Pb_Masks'
spot_folder = 'output'

# Set threshold for total P-Body index, only cells with a total P-Body index larger than 
# The threshold will be used for the filtered results
PB_index_threshold = 1.9



if not path_root:
    raise ValueError('Set path_root before running this script.')
if not conditions:
    raise ValueError('Set at least one condition before running this script.')
if not pbody_chan:
    raise ValueError('Set pbody_chan before running this script.')
if not rna_chan:
    raise ValueError('Set rna_chan before running this script.')

# You only need to run PBody_Charecteristics if you are not running P-Body mRNA metrics.
# The P-Body Charecteristics are contained in the mRNA metrics excel files.
"""
PBody_Charecteristics(
    conditions=conditions,
    diffraction=diffraction,
    path=str(path_root),
    pixel_res=pixel_res,
    pbody_chan=pbody_chan,
    cell_mask_chan=cell_mask_chan,
    nuc_mask_chan=nuc_mask_chan,
    PB_folder=pb_folder,
)
"""

print(f'P-body characteristics workbook: {path_root / f"Pbody Charecteristics {pbody_chan}.xlsx"}')


_, mrna_workbook = run_pbody_mrna_metrics(
    path_root=str(path_root),
    conditions=conditions,
    spot_source=spot_source,
    pbody_dilation_iterations = 2,
    local_shell_iterations = 4,
    pbody_chan=pbody_chan,
    rna_chan=rna_chan,
    cell_mask_chan=cell_mask_chan,
    nuc_mask_chan=nuc_mask_chan,
    segmentation_folder=pb_folder,
    spot_folder=spot_folder,
    diffraction=diffraction,
    pixel_res=pixel_res,
    )

print(f'mRNA metrics workbook: {mrna_workbook}')


# Creates results spreadsheet only containing information from cells that have a total PB Index above
# the P-Body index thresholdp set by the user above, comment out if you do not want filtered results
_, mrna_workbook = run_pbody_mrna_metrics(
    path_root=str(path_root),
    conditions=conditions,
    spot_source=spot_source,
    pbody_dilation_iterations = 2,
    local_shell_iterations = 4,
    pbody_chan=pbody_chan,
    rna_chan=rna_chan,
    cell_mask_chan=cell_mask_chan,
    nuc_mask_chan=nuc_mask_chan,
    segmentation_folder=pb_folder,
    spot_folder=spot_folder,
    diffraction=diffraction,
    pixel_res=pixel_res,
    filter = True,
    PB_index_threshold=PB_index_threshold,
    ) 

print(f'Filtered mRNA metrics workbook: {mrna_workbook}')