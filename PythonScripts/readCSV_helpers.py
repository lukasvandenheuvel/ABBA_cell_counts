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
def plot_starter_cells(brain_df, brain_region_dict, output_path):
    # Starter cell analysis
    starter_cells = brain_df['RAB_TVA']
    starter_cells = starter_cells[starter_cells > 0]
    starter_cells_sorted = sort_hemispheres(starter_cells)
    
    # Plot starter cells
    index = ['%s (%s)'%(brain_region_dict[key],key) for key in starter_cells_sorted.index]
    plt.figure(figsize=(20,5))
    plt.tight_layout()
    b = plt.bar(index, starter_cells_sorted['Sum'])
    t = plt.xticks(rotation=90)
    t = plt.title('Starter cells')
    lbl = plt.ylabel('Starter cells (Rabies+ TVA+)')
    if not(output_path==None):
        output_file = os.path.join(output_path, 'starter_cells.pdf')
        plt.savefig(output_file, bbox_inches='tight')

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
def collect_and_analyze_cell_counts(root, animal_list, tracers, path_to_onotlogy_pickle):
    
    # Store the seperate hemispheres, and the sum of the hemispheres:
    hemispheres = ['Left', 'Right', 'Sum']

    # Load brain ontology (brain hierarchy) --------------------------------------
    with open(path_to_onotlogy_pickle,"rb") as f:
        ontology_dict = pickle.load(f)
    edges = ontology_dict['BrainOntologyEdges']
    tree = ontology_dict['BrainOntologyTree']
    brain_region_dict = ontology_dict['BrainOntologyRegions']

    # Initialize a results dataframe. --------------------------------------------
    # This is a dataframe with hierarchical columns. 
    # Hierarchy: tracer -> animal -> hemisphere.
    iterables = [tracers, animal_list, hemispheres]
    multi_index = pd.MultiIndex.from_product(iterables)
    results = pd.DataFrame(np.nan, index=brain_region_dict.keys(), columns=multi_index)

    # Loop over animals, load the data and normalize counts --------------------
    for animal in animal_list:

        print('Importing slices in '+animal+'...')
        input_path = os.path.join(root, animal, 'results')
        output_path = os.path.join(root, animal, 'results_python')

        # Load regions to exclude for this animal
        path_to_exclusion_file = os.path.join(root, animal, 'RegionsToExclude.csv')
        if not(os.path.exists(path_to_exclusion_file)):
            raise ValueError('Cannot find exclusion file for animal ' + animal + '!')
        exclude_dict = list_regions_to_exclude(path_to_exclusion_file)

        # Load cell counts, excluding the regions we want to exclude
        df_list,slice_regions,slice_data = load_cell_counts(input_path, exclude_dict, edges, tree)
        print('Imported ' + str(len(df_list)) + ' slices.\n')

        # Now comes the tricky part. We'll first concatenate the dataframes
        # of all slices into one big dataframe (brain_df).
        # Then, we combine the rows with the same index (=region name), and sum them.
        # That is, we sum the results (area, cell counts) per region across slices.
        brain_df = pd.concat(df_list)
        brain_df = brain_df.groupby(brain_df.index, axis=0).sum()
        
        # Plot starter cells
        plot_starter_cells(brain_df, brain_region_dict, output_path)

        # Save brain_df
        brain_df.to_csv( os.path.join(output_path, animal+'_cell_counts.csv') )
        print('Raw cell counts are saved to ' + output_path)

        # Normalize the results
        for t in tracers: # loop over tracers ('RAB', 'CTB', ...)

            # Normalize
            normalized_cell_counts = normalize_cell_counts(brain_df, t)

            # Save results per animal
            present_regions = normalized_cell_counts.index.to_list()
            for region in present_regions: # loop over all regions present
                results.loc[region, (t,animal)].update( normalized_cell_counts.loc[region] )

    # Swap hierarchy of columns, to make averaging over animals easier.
    # The new hierarchy will be Tracer -> Hemisphere -> Animal
    results = results.swaplevel(axis=1)
    
    return results

#%%
def average_cell_counts_over_animals(results, tracers):
    # Calculate means and sems -------------------------------------------------
    iterables = [tracers, ['PerHemi', 'SummedHemi'], ['Mean', 'Sem']]
    multi_index = pd.MultiIndex.from_product(iterables)
    mean_results = pd.DataFrame(np.nan, index=results.index, columns=multi_index)

    for t in tracers:

        # Normalization per hemisphere: Treat 'Left' and 'Right' as seperate animals to calculate average
        mean_results.loc[:, (t,'PerHemi','Mean')] = results[t][['Left','Right']].mean(axis=1)
        mean_results.loc[:, (t,'PerHemi','Sem')] = results[t][['Left','Right']].sem(axis=1)

        # Normalization with summed hemispheres: Sum left and right, and average over animals.
        mean_results.loc[:, (t,'SummedHemi','Mean')] = results[t]['Sum'].mean(axis=1)
        mean_results.loc[:, (t,'SummedHemi','Sem')] = results[t]['Sum'].sem(axis=1)
        
    return mean_results
