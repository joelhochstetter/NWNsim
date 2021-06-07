#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This code produces multiple nanowire networks of specified parameters into a folder
# 
# 'multi_generate_networks.py' is designed to be used from terminal (bash, etc.).
# 
#
# Usage:
# python multi_generate_networks.py [argumentName] [argument]
# 
# Saves as a file storing adjacency matrix, wire locations, junction locations, etc. 
#
# Written by Joel Hochstetter


import wires
import numpy as np

import argparse

# Create parser for options
parser = argparse.ArgumentParser(
    description='Handle parameters to generate a network of nanowires and junctions.')


parser.add_argument('--numSims',
    type    = int,
    default = 10,
    help    = 'The number of sims in the network.')


parser.add_argument('--nwires',
    type    = int,
    default = 100,
    help    = 'The number of nanowires in the network.')
    
parser.add_argument('--nwiresMax',
    type    = int,
    default = -1,
    help    = 'The number of nanowires in the network.')


parser.add_argument('--mean_length', 
    type    = float, 
    default = 10.0,
    help    = 'The mean length of the nanowires. Passed to the gamma distribution.')

parser.add_argument('--std_length', 
    type    = float, 
    default = 1.0,
    help    = 'The standard deviation of nanowires length. Passed to the gamma distribution.')

parser.add_argument('--seed',
    type    = int, 
    default = 0,
    help    ='The seed for the random number generator.')

parser.add_argument('--seedMax',
    type    = int, 
    default = -1,
    help    ='The maximum seed for the random number generator.')

parser.add_argument('--Lx',
    type    = float, 
    default = 100,
    help    ='The horizontal length of the network''s physical substrate in micrometres.')

parser.add_argument('--LxMax',
    type    = float, 
    default = -1,
    help    ='The horizontal length of the network''s physical substrate in micrometres.')

parser.add_argument('--Ly',
    type    = float,
    default = -1,
    help    ='The vertical length of the network''s physical substrate in micrometres.')

parser.add_argument('--cent_dispersion', 
    type    = float, 
    default = 700.0,
    help    = 'The width of the generalised normal distribution in units of um.')

parser.add_argument('--shape', 
    type    = float, 
    default = 5.0,
    help    = 'Shape parameter beta. Passed to the generalised normal distribution. Value of 2 is normal. ->inf is uniform.')


parser.add_argument('--density', 
    type    = float, 
    default = -1, #e.g. 0.1
    help    = 'number of wires per um^2.')


parser.add_argument('--folder', 
    type    = str, 
    default = 'connectivity_data',
    help    ='The folder where the output files will be stored.')


parser.add_argument('--oldNameConvention', 
    type    = bool, 
    default = False,
    help    = 'Whether or not to use old name convention (include cent dispersion and date), or new (include lx, ly, no date).')
    
parser.add_argument('--plot', 
    dest    = 'plot_network', 
    action  = 'store_true',
    default = False, 
    help    = 'Flag to plot the figure.')

parser.add_argument('--no-plot', 
    dest    = 'plot_network', 
    action  = 'store_false',
    default = False,
    help    = 'Flag to not plot the figure (default).')    

args = parser.parse_args()

if args.nwiresMax == -1:
    args.nwiresMax = args.nwires
    
if args.seedMax == -1:
    args.seedMax = args.seed + 1

if args.Ly == -1:
    args.Ly = args.Lx
    square = True
else:
    square = False

if args.LxMax == -1:
    args.LxMax = args.Lx

mean_length     = args.mean_length
std_length      = args.std_length
shape           = args.shape
cent_dispersion = args.cent_dispersion
Lx              = args.Lx
Ly              = args.Ly
folder          = args.folder
density         = args.density
oldNameConvention  = args.oldNameConvention

wireList = list(np.unique(np.linspace(args.nwires, args.nwiresMax, args.numSims, dtype = int)))
seedList = range(args.seed, args.seedMax)
LxList   = list(np.unique(np.linspace(args.Lx,     args.LxMax,     args.numSims, dtype = int)))

print('multi_generate_networks:')

if density == -1:
    print('Seeds: ', list(seedList))
    print('Wires: ', list(wireList))
    print('Sizes: ', list(LxList))
    fixedDensity = False
else:
    wireList = [-1]
    print('Density: ', density)
    print('Seeds: ', list(seedList))
    print('Sizes: ', list(LxList))    
    fixedDensity = True

print('Fixed density = ', fixedDensity)
print('Square = ', square)

if oldNameConvention:
    print('Using old name convention')


for nwires in wireList: 
    for seed in seedList:
        for Lx in LxList:
            if square:
                Ly = Lx         
            if fixedDensity:
                nwires = int(density*Lx*Ly)
            
            print('Now generating: nwires = ', nwires, ', seed = ', seed, ' grid = ', Lx, 'x', Ly, 'um^2')
            # Generate the network
            wires_dict = wires.generate_wires_distribution(number_of_wires = nwires,
                                                     wire_av_length = mean_length,
                                                     wire_dispersion = std_length,
                                                     gennorm_shape = shape,
                                                     centroid_dispersion= cent_dispersion,
                                                     this_seed = seed,
                                                     Lx = Lx,
                                                     Ly = Ly, 
                                                     oldNameConvention = oldNameConvention)


            # Get junctions list and their positions
            wires_dict = wires.detect_junctions(wires_dict)

            # Genreate graph object and adjacency matrix
            wires_dict = wires.generate_graph(wires_dict)
            
            
            if not wires.check_connectedness(wires_dict):
                wires_dict = wires.select_largest_component(wires_dict)
            
            
            #Calculate network statistics
            #wires_dict = wires.analyse_network(wires_dict)

            wires.export_to_matlab(wires_dict, folder = folder)
            
            if args.plot_network:
             
                # Plotting tools
                from matplotlib.lines import Line2D
                from matplotlib.patches import Rectangle
                import matplotlib.pyplot as plt

                # Plot pretty pictures of what we just did
                fig, ax = plt.subplots()
                fig.set_size_inches(5,5)

                Lx = wires_dict['length_x']
                Ly = wires_dict['length_y']

                ax.add_patch(Rectangle((0,0), Lx, Ly, color=(1.0, 0.918, 0.0), alpha=0.77))     
                ax = wires.draw_wires(ax, wires_dict)
                ax = wires.draw_junctions(ax, wires_dict)
                ax.set_aspect(1) # set aspect ratio to 1
                ax.set_xlabel(r'x [$\mu$m]')
                ax.set_ylabel(r'y [$\mu$m]')
                ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                ax.axis([-.1*Lx,1.1*Lx,-.1*Lx,1.1*Lx]) # add some space around the unit square
                ax.set_title('Nanowires distribution')
                ax.grid()
                plt.show()       
 
