#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 09:13:52 2021

Initialize a CSV file in which you can indicate which regions to exclude in each slice.

@author: lukasvandenheuvel
"""

import os
import pandas as pd


#%% ------------------------------ SET PARAMETERS ----------------------------
# ============================================================================
path_to_animal = '/Users/lukasvandenheuvel/Documents/GRAFF Lab/2021_RabiesTracing/TRIO/TRIO2_11261_Lukas'

#%% -------------------------------- START SCRIPT ----------------------------
# ============================================================================

output_file = os.path.join(path_to_animal, 'RegionsToExclude.csv')
path_to_csv_files = os.path.join(path_to_animal,'results')

if not(os.path.exists(path_to_csv_files)):
    raise ValueError('No results directory found!')
    
if os.path.exists(output_file):
    ans = input('RegionsToExclude.csv already exists! Are you sure you want to continue and overwrite the existing file? (y/n) ')
    if ans != 'y':
        raise ValueError('User terminated the script.')
    
# Get files
all_files = os.listdir(path_to_csv_files)
# Filter txt files, and remove the '_regions.txt' extension
txt_files = [f for f in all_files if '_regions.txt' in f]

txt_files.sort()

df = pd.DataFrame(txt_files, columns=['Image Name'])
df['Regions to Exclude (Regions may not overlap!)'] = [None] * len(txt_files)
df.to_csv(output_file, index=False)

print('Script finished!')