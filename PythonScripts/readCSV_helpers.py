# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%%
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle

import copy
import json

#%%
def get_image_names_in_folder(path):
    '''
    Returns a list of all files in the directory
    that have a '.txt' extension in them. It removes '_LEFT' and '_RIGHT' from the names.
    '''
    all_files = os.listdir(path)
    # Filter txt files, and remove the '_regions.txt' extension
    txt_files = [f.replace('_regions.txt', '') for f in all_files if '_regions.txt' in f]
    # Get rid of '_LEFT' and '_RIGHT':
    files = []
    for f in txt_files:
        if '_LEFT' in f:
            files.append(f.replace('_LEFT', ''))
        elif '_RIGHT' in f:
            files.append(f.replace('_RIGHT', ''))
        else:
            files.append(f)
            
    # Remove doubles and sort:
    files = list(dict.fromkeys(files))
    files.sort()
    
    return files

#%%
def remove_hemisphere(data, hemisphere):
    '''
    This function removes all regions specific to
    a certain hemisphere (either 'Left' or 'Right')
    from a dataframe.
    
    Inputs
    ------
    data (pandas dataframe)
    A dataframe with the counting results of a slice.
    
    hemisphere (string)
    Choose 'Left' or 'Right'.
    '''
    
    if not(hemisphere=='Left' or hemisphere=='Right'):
        raise ValueError('Hemisphere should be either "Left" or "Right"!')
    
    curr_regs = data.index.tolist()
    for region in curr_regs:
        if hemisphere in region:
            data = data.drop(region, axis=0)
    
    return data

#%% 
def import_txt_file_as_dataframe(path_to_txt, hemisphere):
    '''
    This function reads a txt file into a pandas dataframe.
    It does some additional processing steps to make the handling
    of the data easier in the next steps. These steps are:
    - Create a class name for the Root region called ROOT.
    - Replace NaN values with 0.
    - Convert the Class column to the index of the dataframe.
    '''
    data = pd.read_table(path_to_txt)
    img_name = data.loc[0,'Image Name']
    
    # There is one region (the full slice) called 'Root', which has the class 'NaN'
    # associated to it. Replace this NaN with the more descriptive class root.
    data['Class'] = data['Class'].fillna('wholeroot')

    # Set the region classes as the column index
    data = data.set_index('Class')
    
    # Now remove the 'wholeroot'. We'll use the seperate hemispheres.
    data = data.drop('wholeroot', axis=0)
    
    # If a hemisphere is specified, remove the other hemisphere from the dataframe
    if (hemisphere == 'Left'):
        data = remove_hemisphere(data, 'Right')
    elif (hemisphere == 'Right'):
        data = remove_hemisphere(data, 'Left')
    
    return data,img_name

#%%
def find_region_abbreviation(region_class):
    '''
    This function finds the region abbreviation
    by splitting the class value in the table.
    Example: 'Left: AVA' becomes 'AVA'.
    '''
    try: # try to split the class
        region_abb = region_class.split(': ')[1]
    except: # if splitting gives an error, don't split
        region_abb = str(region_class)
        
    return region_abb

#%%
def filter_uppercase_characters(string):
    '''
    This function returns all uppercase characters in the input string.
    Example: 'ACAd' returns 'ACA'.
    '''
    uppercase = [char for char in string if char.isupper()]
    return ''.join(uppercase)

#%%
def find_regions_and_classes_in_slice(data):
    '''
    This function reads a dataframe that corresponds to a brain slice,
    and returns a dictionary with the names of the classes appearing
    in the slice as keys, and the full region names as the corresponding value.
    Example: 
    region_dict['ACAd'] = 'Anterior cingulate area, dorsal part' 
    '''
    # Get the number of rows, columns in dataframe
    num_regions, num_measurements = data.shape
    # Initialize region_dict
    region_dict = {}
    # Loop through all regions
    for region_class in data.index.tolist():
        # Put the region class name (abbreviation)
        # as key in the dictionary, and the full region name as corresponding value.
        region_name = data.loc[region_class,'Name']
        region_class = find_region_abbreviation(region_class)
        region_dict[region_class] = region_name
        
    return region_dict

#%%
def list_all_subregions(region_to_list, tree):
    '''
    This function lists all subregions belonging to region_to_list.
    It includes subregions on all hierarchical levels.
    
    Inputs
    ------
        region_to_list (str)
        Name of the region of which you want to know the subregions.
        
        tree (dict)
        Dictionary with brain regions as keys and their children as value.
        
    Output
    ------
        subregions (list)
        List of all subregions.
    '''
    
    d = copy.deepcopy(tree) 
    subregions = [region_to_list]

    # Initialize the region to list as child.
    child = region_to_list

    # Keep walking as long as the region_to_list is in the
    # copied brain dictionary.
    while region_to_list in d.keys(): 

        region = child
        # Add this subregion to the list, if it wasn't in there already
        if region not in subregions:
            subregions.append(region)

        if region in d.keys():   # If the region still has unlisted children ...
            parent = region      # the region now becomes a parent ...
            child = d[parent][0] # and we take the first of its children as new child.

        else:                         # If the region has no unlisted children ...
            d[parent].remove(region)  # we remove it from its parent ...
            child = region_to_list    # and go back to start.
            if len(d[parent])==0:     # We remove the parent if this was the last child.
                del d[parent]
    
    return subregions

#%%
def list_regions_to_exclude(path_to_exclusion_file):
    '''
    Read the csv file containing the regions to exclude for each image,
    and summarize the information in a dictionary.
    The exclusions_file is to be initialized with initExclusionFile
    '''
    exclude_df = pd.read_csv(path_to_exclusion_file, sep=r'[;,]', index_col='Image Name', engine='python')
    exclude_dict = {}

    for img in exclude_df.index:
        
        exclude_dict[img] = []
        to_exclude = exclude_df['Regions to Exclude (Regions may not overlap!)'].loc[img]
        if type(to_exclude) == str: # if there are regions to exclude
            for region in to_exclude.split('/ '):
                
                # If no hemisphere was specified, add both hemispheres:
                if not 'Right' in region and not 'Left' in region:
                    exclude_dict[img].append('Right: '+region)
                    exclude_dict[img].append('Left: '+region)
                # If left or right hemisphere was specified:
                else:
                    exclude_dict[img].append(region)
                
    return exclude_dict
    
#%%
def exclude_regions(df, regs_to_exclude, edges, tree):
    '''
    Take care of regions to be excluded from the analysis.
    If a region is to be excluded, 2 things must happen:
    (1) The cell counts of that region must be subtracted from all
        its parent regions.
    (2) The region must disappear from the data, together with all 
         its daughter regions.
    '''

    # Take care of regions to exclude
    for reg_hemi in regs_to_exclude:
        hemi = reg_hemi.split(': ')[0]
        reg = reg_hemi.split(': ')[1]

        # Step 1: subtract counting results of the regions to be excluded
        # from their parent regions.
        child = reg
        while True:
            parent = edges[child]
            row = hemi+': '+parent
            # Subtract the counting results from the parent region.
            # Use fill_value=0 to prevent "3-NaN=NaN".
            df.loc[row] = df.loc[row].subtract( df.loc[reg_hemi], fill_value=0 )
            if (parent == 'root'):
                break
            child = parent # go one step further in the tree

        # Step 2: Remove the regions that should be excluded
        # together with their daughter regions.
        subregions = list_all_subregions(reg, tree)
        for subreg in subregions:
            row = hemi+': '+subreg
            if row in df.index:
                df = df.drop(row)
            
    return df

#%%
def load_cell_counts(root, exclude_dict, edges, tree):
    '''
    Function to load cell counts, stored in .csv files in the 'root' directory,
    as Pandas dataframes.
    '''
    
    # Get the image names present in root (e.g. "Image_01.vsi - 10x_01")
    # and the names of all files present in root (e.g. "Image_01.vsi - 10x_01 LEFT_regions.txt")
    img_names = get_image_names_in_folder(root)
    file_names = os.listdir(root)
    
    # Init dicts and lists to store data
    slice_regions = {}   # which regions do we have per slice?
    slice_data = {}      # what are the cell counts per slice?
    df_list = []         # list of all slice dataframes
    
    # Loop through the image names
    for f in img_names:

        # The following variables will be used to find out whether we have seperate files
        # for seperate hemispheres, or just one file containing both hemispheres.
        fname = f + '_regions.txt'
        fname_left = f + '_LEFT' + '_regions.txt'
        fname_right = f + '_RIGHT' + '_regions.txt'
        both_hemi = False
        right_hemi = False
        left_hemi = False
        regs_to_exclude = []

        # Read text file into a Pandas dataframe
        if fname_left in file_names: # if we have img_name LEFT_regions.txt in folder
            left_hemi = True
            path = os.path.join(root, fname_left)
            data_left,img_name_left = import_txt_file_as_dataframe(path, 'Left')
            regs_to_exclude = regs_to_exclude + exclude_dict[fname_left]
        if fname_right in file_names: # if we have img_name RIGHT_regions.txt in folder
            right_hemi = True
            path = os.path.join(root, fname_right)
            data_right,img_name_right = import_txt_file_as_dataframe(path, 'Right')
            regs_to_exclude = regs_to_exclude + exclude_dict[fname_right]
        if fname in file_names:       # if we have img_name_regions.txt (no hemisphere specification)
            both_hemi = True
            path = os.path.join(root, fname)
            data,img_name = import_txt_file_as_dataframe(path, 'Both')
            regs_to_exclude = regs_to_exclude + exclude_dict[fname]

        # Check for safety: we either have ONE file for both hemispheres,
        # or (max 2) file(s) for seperate hemispheres. Else, raise and error.
        if (left_hemi and both_hemi) or (right_hemi and both_hemi):
            raise ValueError('Either LEFT and/or RIGHT, or no hemisphere specification. But not both!')
        if not(left_hemi) and not(right_hemi) and not(both_hemi):
            raise ValueError('Filename not found!')

        # Combine left and right, if they were both present
        if left_hemi and right_hemi:        # if we have both left and right, combine dataframes
            data = pd.concat([data_left, data_right])
        elif left_hemi and not(right_hemi): # if we have only left, data = data_left
            data = data_left
        elif not(left_hemi) and right_hemi: # if we have only right, data = data_right
            data = data_right

        # Find regions in current slice
        region_dict = find_regions_and_classes_in_slice(data)

        # Combine cell counts
        df = sum_cell_counts(data)
        
        # Take care of regions to be excluded
        df = exclude_regions(df, regs_to_exclude, edges, tree)
        
        # Store results in dictionaries / lists
        slice_regions[f] = region_dict
        slice_data[f] = df
        df_list.append(df)
    
    return df_list,slice_regions,slice_data
    
#%%
def sum_cell_counts(data):
    '''
    This function takes as input raw data from a csv file (data = a dataframe created with pd.read_csv).
    It converts counts (Num CTB (only), Num CTB: Rabies (only), etc) to number of detected cells ('CTB, RAB', etc).
    To do this it sums the relevant counts.
    '''
    
    # The parameters we are interested in
    param_list = ['area',                             # DAPI area
                  'CTB', 'RAB', 'TVA',                # single positives + double positives + triple positives
                  'CTB_RAB', 'CTB_TVA', 'RAB_TVA',    # double positives + triple positives
                  'CTB_RAB_TVA']                      # triple positives
    
    # Make an empty table with rows = all regions in current slice, columns = param_list
    df = pd.DataFrame(np.nan, index=data.index, columns=param_list)
    
    # Warning: below, you'll notice that columns are summed a bit weirdly.
    # I used df['c'] = df[['a','b']].sum(axis=1, min_count=1) to sum up columns 'a' and 'b'.
    # min_count=1 ensures that the sum of NaN values is NaN (and not 0).

    df['area'] = data['DAPI: DAPI area um^2']
    # single positives + double positives + triple positives
    df['CTB'] = data[['Num CTB', 'Num CTB: Rabies', 'Num CTB: TVA', 'Num CTB: Rabies: TVA']].sum(axis=1, min_count=1)
    df['RAB'] = data[['Num Rabies', 'Num CTB: Rabies', 'Num Rabies: TVA', 'Num CTB: Rabies: TVA']].sum(axis=1, min_count=1)
    df['TVA'] = data[['Num TVA', 'Num CTB: TVA', 'Num Rabies: TVA', 'Num CTB: Rabies: TVA']].sum(axis=1, min_count=1)
    # double positives + triple positives
    df['CTB_RAB'] = data[['Num CTB: Rabies', 'Num CTB: Rabies: TVA']].sum(axis=1, min_count=1)
    df['CTB_TVA'] = data[['Num CTB: TVA', 'Num CTB: Rabies: TVA']].sum(axis=1, min_count=1)
    df['RAB_TVA'] = data[['Num Rabies: TVA', 'Num CTB: Rabies: TVA']].sum(axis=1, min_count=1)
    # triple positives
    df['CTB_RAB_TVA'] = data['Num CTB: Rabies: TVA']  

    # Return only those regions where DAPI was found
    return df[df['area'] > 0]

#%%
def init_dict(key_list, init_value):
    '''
    This function initializes a dictionary with keys as specified in key_list,
    and corresponding values in init_value.
    
    key_list: list with keys
    init_value: value to initialize with
    '''
    output_dict = {}
    for key in key_list:
        output_dict[key] = init_value
    
    return output_dict

#%%
def sort_hemispheres(data):
    '''
    Function takes as input a dataframe with only one column (the data to plot).
    The rows represent regions, with left and right seperated.
    
    The output is a dataframe with three columns:
    'Left' for each region, 'Right' for each region and the sum of the two.
    '''
    
    all_regions = data.index.to_list()
    present_regions = [find_region_abbreviation(r) for r in all_regions]
    present_regions = list(dict.fromkeys(present_regions)) # remove doubles
    
    # Put normalized counts of seperate hemispheres in seperate columns
    data_sorted = pd.DataFrame(np.nan, index=present_regions, columns=['Left', 'Right', 'Sum'])
    for region in present_regions:
        if 'Left: ' + region in all_regions:
            data_sorted['Left'][region] = data['Left: ' + region]
        if 'Right: ' + region in all_regions:
            data_sorted['Right'][region] = data['Right: ' + region]

    # Sum of left and right
    data_sorted['Sum'] = data_sorted[['Left','Right']].sum(axis=1, min_count=1)
    
    return data_sorted

#%%
def normalize_cell_counts(brain_df, tracer):
    '''
    Do normalization of the cell counts for one tracer.
    The tracer can be any column name of brain_df, e.g. 'TVA', 'RAB', 'CTB_RAB'.
    The output will be a dataframe with three columns: 'Left', 'Right' and 'Sum'.
    The columns 'Left' and 'Right' are normalized w.r.t one hemisphere only.
    The 'Sum' column is normalized w.r.t the whole brain.
    '''
    
    # Sort the hemispheres (Get 'Left', 'Right' and 'sum' as seperate columns)
    area = sort_hemispheres(brain_df['area'])
    cell_counts = sort_hemispheres(brain_df[tracer])

    # Get the the brainwide area and cell counts (corresponding to the root)
    brainwide_area = area.loc['root']
    brainwide_cell_counts = cell_counts.loc['root']

    # Do the normalization for each column seperately.
    norm_cell_counts = (cell_counts / area) / (brainwide_cell_counts / brainwide_area)
    
    return norm_cell_counts
    
#%%
def plot_bidirectional_bar_chart(data, x_label, brain_region_dict, errorbars=None):
    '''
    Plot results from different hemispheres as bidirectional bar chart.
    '''
    
    if not(data.shape[1]==3):
        raise ValueError('The dataframe to plot should have 3 columns, one for left and one for right and one for sum.')
        
    font_color = '#525252'
    hfont = {'fontname':'Calibri'}
    fontsize = 35
    facecolor = '#eaeaf2'
    color_red = '#fd625e'
    color_blue = '#01b8aa'
    index = ['%s (%s)'%(brain_region_dict[key],key) for key in data.index]
    column_left = data['Left']
    column_right = data['Right']
    xerr_left = None if errorbars is None else errorbars['Left']
    xerr_right = None if errorbars is None else errorbars['Right']
    title_left = 'Left'
    title_right = 'Right'
    max_value  = np.nanmax(data.to_numpy())
    
    # Generate subplots
    fig, axes = plt.subplots(figsize=(50,120), facecolor=facecolor, ncols=2, sharey=True)
    fig.tight_layout()
    
    # Plot bars and vertical line at x=1
    axes[0].barh(index, column_left, align='center', color=color_red, zorder=10, xerr=xerr_left)
    axes[0].set_title(title_left, fontsize=fontsize, pad=15, color=color_red, **hfont)
    axes[0].axvline(x=1, c='k', linestyle='--')
    axes[0].set_xlim([0, max_value])
    axes[1].barh(index, column_right, align='center', color=color_blue, zorder=10, xerr=xerr_right)
    axes[1].set_title(title_right, fontsize=fontsize, pad=15, color=color_blue, **hfont)
    axes[1].set_xlim([0, max_value])
    axes[1].axvline(x=1, c='k', linestyle='--')

    # Axis specifications
    axes[0].invert_xaxis() 
    axes[1].set_xticks(axes[0].get_xticks())

    # Adjust subplots to fit them next to each other and add xlabel
    plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
    xlabel = fig.text(0.565, 0.09, x_label, ha='center', fontsize=fontsize)
    
    return fig

#%%
def plot_results(data, x_label):
    
    # Collect nonzero counts (for plotting) -------------------------
    nonzero_data = data[data > 0]
    
    # Sort by sum of left and right
    data_sorted = sort_hemispheres(nonzero_data)
    data_to_plot = data_sorted.sort_values(by=['Sum'])
    data_to_plot = data_to_plot.drop('Sum', axis=1) # get rid of 'Sum' column, it was only used for sorting

    # Plot results -----------------------------------
    fig = plot_bidirectional_bar_chart(data_to_plot, x_label)
    
    return fig

#%%
def plot_horizontal_bar_chart(data, brain_region_dict):
    
    # remove all regions without any cells present, and sort by sum.
    data = data[data['Mean'] > 0].sort_values(by=['Mean'])

    mean = data['Mean']
    errorbars = data['Sem']

    font_color = '#525252'
    fontsize = 35
    facecolor = '#eaeaf2'
    color_red = '#fd625e'
    color_blue = '#01b8aa'
    index = ['%s (%s)'%(brain_region_dict[key],key) for key in data.index]
    xerr = None if errorbars is None else errorbars
    max_value  = np.nanmax(mean.to_numpy()) + np.nanmax(errorbars.to_numpy())

    # Generate figure
    fig = plt.figure(figsize=(50,120), facecolor=facecolor)
    fig.tight_layout()

    # Plot bars and vertical line at x=1
    plt.barh(index, mean, align='center', color=color_red, zorder=10, xerr=xerr)
    plt.axvline(x=1, c='k', linestyle='--')
    plt.xlim([0, 1.2 * max_value])
    
    return fig
