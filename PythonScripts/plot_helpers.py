#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 17:15:20 2022

@author: lukasvandenheuvel
"""

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import networkx as nx
import numpy as np
import os

from readCSV_helpers import sort_hemispheres

def plot_plotly_graph(G,pos):
    '''
    This function plots the brain hierarchy using plotly.
    G = networkx graph
    pos = node positions (as a dictionary)
    '''

    nx.set_node_attributes(G, pos, 'pos')

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            color=[],
            size=5,
            line_width=.1))
    
    node_colors = []
    node_text = []
    for node_id in G.nodes():
        # Number of connections as color
        node_colors.append('#'+G.nodes()[node_id]['color_hex_triplet'])
        # Region name as text to show
        node_text.append(G.nodes()[node_id]['region_name'] + '(' +
                         G.nodes()[node_id]['acronym'] + ')') 

    node_trace.marker.color = node_colors
    node_trace.text = node_text
    
    fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='Hierarchy of brain regions',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
    return fig

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
