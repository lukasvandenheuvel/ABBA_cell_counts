#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 08:31:08 2021

Statistical quantification of tracing in- and outputs.
Outputs results in a folder called 'results_python' that is created inside the root.

@author: lukasvandenheuvel
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import copy
import pickle

from readCSV_helpers import *

#%% ------------------------------ SET PARAMETERS ----------------------------
# ============================================================================

animal_list = ['TRIO2_11876_Lukas_v2']
root = '/Users/lukasvandenheuvel/Documents/GRAFF Lab/2021_RabiesTracing/TRIO/'
path_to_onotlogy_pickle = '../AllenMouseBrainOntology.pk'

tracers = ['RAB', 'CTB', 'TVA']         # Tracers we are interested in.

#%% -------------------------------- START SCRIPT ----------------------------
# ============================================================================

hemispheres = ['Left', 'Right', 'Sum']  # Store the seperate hemispheres, and the sum of the hemispheres.

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

#%% Loop over animals, load the data and normalize counts --------------------
for animal in animal_list:

    print('Importing slices in '+animal+'...')
    input_path = os.path.join(root, animal, 'results')

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
    starter_cells = brain_df['RAB_TVA']
    starter_cells = starter_cells[starter_cells > 0]
    starter_cells_sorted = sort_hemispheres(starter_cells)
    
    # Plot starter cells and save brain_df
    index = ['%s (%s)'%(brain_region_dict[key],key) for key in starter_cells_sorted.index]
    plt.figure(figsize=(20,5))
    plt.tight_layout()
    b = plt.bar(index, starter_cells_sorted['Sum'])
    t = plt.xticks(rotation=90)
    t = plt.title('Starter cells')
    lbl = plt.ylabel('Starter cells (Rabies+ TVA+)')
    output_path = os.path.join(root, animal, 'results_python')
    if not(os.path.isdir(output_path)):
        os.mkdir(output_path)
    output_file = os.path.join(output_path, animal+'_starter_cells.pdf')
    plt.savefig(output_file, bbox_inches='tight')
    brain_df.to_csv( os.path.join(output_path, animal+'_cell_counts.csv') )

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

#%% Calculate means and sems -------------------------------------------------
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

#%% Save and plot results -----------------------------------------------------
output_path = os.path.join(root, 'results_python')
if not(os.path.exists(output_path)):
    os.mkdir(output_path)
    print('\nCreated a new results_python folder in ' + root + '\n')
else:
    print('\n! A results_python folder already existed in root. I am overwriting previous results!\n')

results.to_csv( os.path.join(output_path, 'results_cell_counts.csv') )
mean_results.to_csv( os.path.join(output_path, 'results_mean_cell_counts.csv') )
print('\nGenerating plots ...')

#%% Plot Rabies+ normalized per hemisphere


x_label = '(Rabies+ / Dapi area) / (total Rabies+ / total Dapi area)'

rabies = mean_results['RAB']['PerHemi']
fig = plot_horizontal_bar_chart(rabies, brain_region_dict)

plt.title('Rabies+ normalized per hemisphere',fontsize=35)
output_file = os.path.join(output_path, 'rabies_normalized_per_hemi.pdf')
plt.savefig(output_file, bbox_inches='tight')
plt.show()

#%% Plot Rabies+ normalized with summed hemispheres
rabies = mean_results['RAB']['SummedHemi']
fig = plot_horizontal_bar_chart(rabies, brain_region_dict)

plt.title('Rabies+ normalized as summed hemispheres',fontsize=35)
output_file = os.path.join(output_path, 'rabies_normalized_sum.pdf')
plt.savefig(output_file, bbox_inches='tight')
plt.show()

#%% Plot TVA+ normalized per hemisphere
x_label = '(TVA+ / Dapi area) / (total TVA+ / total Dapi area)'

tva = mean_results['TVA']['PerHemi']
fig = plot_horizontal_bar_chart(tva, brain_region_dict)

plt.title('TVA+ normalized per hemisphere',fontsize=35)
output_file = os.path.join(output_path, 'tva_normalized_per_hemi.pdf')
plt.savefig(output_file, bbox_inches='tight')
plt.show()

#%% Plot TVA+ normalized with summed hemispheres
tva = mean_results['TVA']['SummedHemi']
fig = plot_horizontal_bar_chart(tva, brain_region_dict)

plt.title('TVA+ normalized as summed hemispheres',fontsize=35)
output_file = os.path.join(output_path, 'tva_normalized_sum.pdf')
plt.savefig(output_file, bbox_inches='tight')
plt.show()

#%% Plot CTB+ normalized per hemisphere
ctb = mean_results['CTB']['PerHemi']
fig = plot_horizontal_bar_chart(ctb, brain_region_dict)

plt.title('CTB+ normalized per hemisphere',fontsize=35)
output_file = os.path.join(output_path, 'ctb_normalized_per_hemi.pdf')
plt.savefig(output_file, bbox_inches='tight')
plt.show()

#%% Plot CTB+ normalized with summed hemispheres
ctb = mean_results['CTB']['SummedHemi']
fig = plot_horizontal_bar_chart(ctb, brain_region_dict)

plt.title('CTB+ normalized as summed hemispheres',fontsize=35)
output_file = os.path.join(output_path, 'ctb_normalized_sum.pdf')
plt.savefig(output_file, bbox_inches='tight')
plt.show()

print('Script finished!')
