#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module generates a distrbution of nanowires on 2D domain, akin to the 
where atomic switch networks are grown. 

The basic process consists in choosing a random center point for the wire in 
the unit square and then chooses a random angle \theta \in (0,\pi) as the 
wire's orientation.

.. moduleauthor:: Paula Sanz-Leon
.. moduleauthor:: Miro Astore
"""

from __future__ import  division

from itertools import *
from scipy.io import savemat
from scipy.spatial.distance import cdist
from scipy.stats import gennorm

import numpy as np
import networkx as nx
import time
import pickle

from scipy.stats import gennorm

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

def generate_wires_distribution(number_of_wires=1500, 
                                wire_av_length=14.0, 
                                wire_dispersion=5.0,
                                centroid_dispersion=1200.0,  
                                gennorm_shape = 5,
                                Lx=3e3, Ly=3e3, 
                                this_seed=42,
                                oldNameConvention = False):
    '''
    Drops nanowires on the device of sides Lx, Ly. 
    
    Parameters
    ----------
    number_of_wires : int 
        Total number of wires to be sampled
    wire_av_length : float 
        Average wire length in mum
    wire_dispersion : float 
        Dispersion/scale of length distribution in mum
    wire_length : float 
        Length of the nanowire in mum (default = 14)
    centroid_dispersion : float 
        Scale parameter for the general normal distribution from 
        which centroids of wires are drawn in mum
    gennorm_shape : float 
        Shape parameter of the general normal distribution from 
        which centroids of wires are drawn. As this number increases, 
        the distribution approximates a uniform distribution.
    Lx : float 
        Horizontal legth of the device in mum
    Ly : float 
        Vertical length of the device in mum
    seed : int
        Seed of the random number generator to always generate the same distribution
    
    Returns
    -------
    dict
        A dictionary with the centre coordinates, the end point coordinates, and
        orientations. The `outside` key in the dictionary is 1 when
        the wire intersects an edge of the device and is 0 otherwise.

    '''

    np.random.seed(this_seed)
    
    # wire lengths
    wire_lengths     = generate_dist_lengths(number_of_wires, wire_av_length, wire_dispersion)
#    xc, yc = generate_dist_centroids(number_of_wires, loc=Lx/2, scale=centroid_dispersion, beta=gennorm_shape)
    xc, yc = np.random.rand(number_of_wires) * Lx, np.random.rand(number_of_wires) * Ly
    theta  = generate_dist_orientations(number_of_wires)
    
    xa, ya  = xc - wire_lengths/2.0 * np.cos(theta), yc - wire_lengths/2.0 * np.sin(theta) # coordinates for one end 
    xb, yb  = xc + wire_lengths/2.0 * np.cos(theta), yc + wire_lengths/2.0 * np.sin(theta) # coordinates for the other end

    # Compute the Euclidean distance between pair of region centres
    wire_distances = cdist(np.array([xc, yc]).T, np.array([xc, yc]).T, metric='euclidean')    

    # Find values outside the domain
    a = np.where(np.vstack([xa, xb, ya, yb]) < 0.0, True, False).sum(axis=0)
    b = np.where(np.vstack([xa, xb]) > Lx, True, False).sum(axis=0)
    c = np.where(np.vstack([ya, yb]) > Ly, True, False).sum(axis=0)
    
    outside = a + b + c

    return dict(xa=xa,ya=ya,
                xc=xc,yc=yc,
                xb=xb,yb=yb,
                theta=theta,
                avg_length = wire_av_length,
                dispersion = wire_dispersion,
                centroid_dispersion = centroid_dispersion,
                gennorm_shape = gennorm_shape,
                this_seed = this_seed,
                outside = outside,
                length_x = Lx,
                length_y = Ly,
                number_of_wires = number_of_wires,
                wire_distances = wire_distances,
                oldNameConvention = oldNameConvention)


def generate_dist_centroids(number_of_wires, loc = 10, scale = 50, beta=5):
    '''
    Generates the 2D coordinates from a general normal distribution
    '''

    # This is a uniform random distribution
    xc = np.random.rand(number_of_wires) * Lx # uniform random in [0,Lx) for each coordinate dimension
    yc = np.random.rand(number_of_wires) * Ly # uniform random in [0,Lx) for each coordinate dimension

    
    #xc = gennorm.rvs(beta, loc=loc, scale=scale, size=int(number_of_wires))
    #yc = gennorm.rvs(beta, loc=loc, scale=scale, size=int(number_of_wires))
    
    return xc, yc


def generate_dist_lengths(number_of_wires, wire_av_length, wire_dispersion):
    '''
    Generates the distribution of wire lengths
    '''

    gamma_shape = (wire_av_length / wire_dispersion)**2
    gamma_scale = wire_dispersion**2 / wire_av_length

    wire_lengths = np.random.gamma(gamma_shape, gamma_scale, int(number_of_wires))

    return wire_lengths


def generate_dist_orientations(number_of_wires):

     return np.random.rand(int(number_of_wires))*np.pi # uniform random angle in [0,pi)


def generate_lattice (number_of_wires=1500, 
                          wire_length=10.0):
    '''
        generates a lattice something like this, approximates a square shape as best as possible given number of wires
        
              |           |
              |           |     
    wire1}----o----{o}----o----{wire 2 
              |           |         o = junction
              |           |         { = end of wire
            wire4        wire3
        
    '''
    # make sure number of wires is even
    if number_of_wires % 2 != 0:
        number_of_wires = number_of_wires+1

#   find squarest way to arrange grid, arrange side lengths such that second number is higher, this will be the x side length
    ideal_square = int(np.sqrt(number_of_wires/2))
    factors = [0,0]
    for i in range(ideal_square,int(number_of_wires/2)+1):
        if (number_of_wires/2) % i == 0:
            factors[0] = int((number_of_wires/2)/i)
            factors[1] = i
            break

    factors = sorted(factors)
    
    #arrange points in grid
    float_range_x = range(factors[1])
    float_range_y = range(factors[0])

    float_range_x = [float(x) for x in float_range_x]
    float_range_y = [float(x) for x in float_range_y]

    float_range_x = [x+0.5 for x in float_range_x]
    float_range_y = [y+0.5 for y in float_range_y]

    xc = [i*wire_length for i in float_range_x]*factors[0]*2
    yc = [i*wire_length for i in float_range_y]*factors[1]*2
    
    xc = np.asarray(xc)
    yc = np.asarray(yc) 
    theta = np.zeros(int(number_of_wires/2))
    temparr = np.ones(int(number_of_wires/2))*np.pi/2
    theta = np.append(theta,temparr)
    
    xa, ya = xc - wire_length/(2.0) * np.cos(theta), yc - wire_length/(2.0) * np.sin(theta) # coordinates for one end 
    xb, yb = xc + wire_length/(2.0) * np.cos(theta), yc + wire_length/(2.0) * np.sin(theta) # coordinates for the other end
    
    Lx = max(np.append(xa,xb))
    Ly = max(np.append(ya,yb))
    
    wire_distances = cdist(np.array([xc, yc]).T, np.array([xc, yc]).T, metric='euclidean')

    # Find values outside the domain
    a = np.where(np.vstack([xa, xb, ya, yb]) < 0.0, True, False).sum(axis=0)
    b = np.where(np.vstack([xa, xb]) > Lx, True, False).sum(axis=0)
    c = np.where(np.vstack([ya, yb]) > Ly, True, False).sum(axis=0)

    outside = a + b + c
    
    return      dict(xa=xa,ya=ya,
                xc=xc,yc=yc,
                xb=xb,yb=yb,
                theta=theta,
                avg_length = wire_length,
                centroid_dispersion = 0,
                this_seed = 0,
                outside = outside,
                length_x = Lx,
                length_y = Ly,
                number_of_wires = number_of_wires,
                wire_distances = wire_distances,
                gennorm_shape = 0,
                dispersion = 0,
                factors = factors) # What is inside factors?


def find_segment_intersection(p0, p1, p2, p3):
    """
    Find *line segments* intersection using line equations and 
    some boundary conditions.

    First segment is defined between p0, p1 and 
    second segment is defined between p2, p3
          p2
          |  
    p0 ------- p1
          |
          p3
    Parameters
    ----------
    p0 : array
        x, y coordinates of first wire's start point 
    p1 : array
        x, y coordinates of first wire's end point
    p2 : array
        x, y coordinates of second wire's start point 
    p3 : array
        x, y coordinates of second wire's end point

    Returns
    -------
    xi, yi: float 
       x, y coordinates of the intersection

    TODO: + change input to a list instead of individual points; or,
          + make point a class with x, y coordinates so we avoid using 
          indexing (x: pX[0]; y:pX[1])
          + polish these docstring with standard input/ouput definitions
    """
    
    # Check that points are not the same
    if np.array_equal(p0, p1) or np.array_equal(p2, p3):
        return False 

     # Check that an overlapping interval exists
    if max(p0[0], p1[0]) < min(p2[0], p3[0]) or max(p2[0], p3[0]) < min(p0[0], p1[0]):
        return False 
    else:
        # xi, yi have to be included in
        interval_xi = [max(min(p0[0],p1[0]), min(p2[0],p3[0])), min(max(p0[0],p1[0]), max(p2[0],p3[0]))]
        interval_yi = [max(min(p0[1],p1[1]), min(p2[1],p3[1])), min(max(p0[1],p1[1]), max(p2[1],p3[1]))]

    # Find the intersection point between nanowires
    A1 = (p0[1]-p1[1])/(p0[0]-p1[0]) # will fail if division by zero
    A2 = (p2[1]-p3[1])/(p2[0]-p3[0]) 
    b1 = p0[1] - A1*p0[0] 
    b2 = p2[1] - A2*p2[0] 

    xi = (b2 - b1) / (A1 - A2)
    yi = A1 * xi + b1

    #The last thing to do is check that xi, yi are included in interval_i:
    if xi  < min(interval_xi) or xi > max(interval_xi):
        return False 
    elif yi < min(interval_yi) or yi > max(interval_yi):
        return False
    else:
        return xi, yi


def select_largest_component(wires_dict):
    """
    Find and select largest connected component of the original graph G.
    Throws away unconnected components and updates all the keys in wires_dict 
    """
    
    wires_dict['G'] = max(nx.connected_component_subgraphs(wires_dict['G']), key=len)
    nw = len(wires_dict['G'].nodes())
    nj = len(wires_dict['G'].edges())   
    
    logging.info("The largest component has %5d nodes and %6d edges", nw, nj)

    # Replace values in the dictionary
    wires_dict['generating_number_of_wires']     =  wires_dict['number_of_wires']    
    wires_dict['number_of_wires']     = nw
    wires_dict['number_of_junctions'] = nj
    wires_dict['xa']    = wires_dict['xa'][sorted(wires_dict['G'].nodes())] 
    wires_dict['ya']    = wires_dict['ya'][sorted(wires_dict['G'].nodes())] 
    wires_dict['xb']    = wires_dict['xb'][sorted(wires_dict['G'].nodes())] 
    wires_dict['yb']    = wires_dict['yb'][sorted(wires_dict['G'].nodes())]
    wires_dict['xc']    = wires_dict['xc'][sorted(wires_dict['G'].nodes())] 
    wires_dict['yc']    = wires_dict['yc'][sorted(wires_dict['G'].nodes())]
    wires_dict['theta'] = wires_dict['theta'][sorted(wires_dict['G'].nodes())]
        			
    # Keep old edge_list
    old_edge_list = [(ii, kk) for ii, kk in  zip(wires_dict['edge_list'][:, 0], wires_dict['edge_list'][:, 1])]
    # Remove old edge list
    wires_dict = remove_key(wires_dict, 'edge_list') 
    # Save indices of intersections in the old graph
    ind_dict = {key:value for value,key in enumerate(old_edge_list)}
    new_edge_list = sorted([kk if kk[0] < kk[1] else (kk[1], kk[0]) for kk in wires_dict['G'].edges()], key=lambda x: x[0])
    # Find intersection between the two sets
    inter = set(ind_dict).intersection(new_edge_list)
    # Retrieve edge indices/positions from the old list
    edges_idx = [ind_dict[idx] for idx in inter]
       
    # These have length equal to number of junctions -- only get the ones we need
    wires_dict['xi'] = wires_dict['xi'][edges_idx] 
    wires_dict['yi'] = wires_dict['yi'][edges_idx] 
    
    # Get contiguous numbering of nodes
    # Build node mapping 
    node_mapping    = {key:value for value, key in enumerate(sorted(wires_dict['G'].nodes()))}
    # This  step also renames edges list
    wires_dict['G'] =  nx.relabel_nodes(wires_dict['G'] , node_mapping)

    # Swap node vertices if vertex 0 is larger than vertex 1, then sort by first element
    wires_dict['edge_list'] = np.asarray(sorted([kk if kk[0] < kk[1] else (kk[1], kk[0]) for kk in wires_dict['G'].edges()], key=lambda x: x[0]))
    
    # Save adjacency matrix of new graph
    wires_dict = remove_key(wires_dict, 'adj_matrix') 
    wires_dict = generate_adj_matrix(wires_dict)

    wire_distances = cdist(np.array([wires_dict['xc'], wires_dict['yc']]).T, np.array([wires_dict['xc'], wires_dict['yc']]).T, metric='euclidean')    
    wires_dict['wire_distances'] = wire_distances

    return wires_dict 


def detect_lattice_junctions(wires_dict):

    factors = wires_dict['factors']
    number_of_wires = wires_dict['number_of_wires']
    
    xi, yi, edge_list = [], [], []
    
    #place junctions at the 'cross' between wires'
    for i in range (int(wires_dict['number_of_wires']/2)):
        xi.append(wires_dict['xc'][i])
        yi.append(wires_dict['yc'][i])
        edge_list.append([i, int(i+wires_dict['number_of_wires']/2)])
    
    #Placing junctions along the horizontal
    for_loop_sequence = range(int(number_of_wires/2)-1)
    deleted_indices = range(factors[1]-1,int(number_of_wires/2)-1,factors[1])
    #Skip the last wire on every row
    for i in deleted_indices[:]:
        if i in for_loop_sequence:
            for_loop_sequence.remove(i)
            deleted_indices.remove(i)
    
    for i in for_loop_sequence:
            edge_list.append([i,i+1])
            # place junction coordinates at end of wire
            xi.append(wires_dict['xa'][i])
            yi.append(wires_dict['ya'][i])
            
    # Place junction at points of vertical lattice wires

    for_loop_sequence = range(int(number_of_wires/2), int(number_of_wires)-factors[1])

    for i in for_loop_sequence:
            edge_list.append([i,i+factors[1]])
            xi.append(wires_dict['xb'][i])
            yi.append(wires_dict['yb'][i])

    wires_dict['number_of_junctions'] = len(np.asarray(edge_list))
    wires_dict['xi'] = np.asarray(xi)
    wires_dict['yi'] = np.asarray(yi)
    wires_dict['edge_list'] = np.asarray(edge_list)
    return wires_dict


def detect_junctions(wires_dict):
    """
    Find all the pairwise intersections of the wires contained in wires_dict.
    Adds four keys to the dictionary: junction coordinates, edge list, and
    number of junctions.

    Parameters
    ----------
    wires_dict: dict

    Returns
    -------
    wires_dict: dict 
        with added keys

    """
    logging.info('Detecting junctions')
    xi, yi, edge_list = [], [], []
    for this_wire, that_wire in combinations(range(wires_dict['number_of_wires']), 2):

        xa, ya = wires_dict['xa'][this_wire], wires_dict['ya'][this_wire]
        xb, yb = wires_dict['xb'][this_wire], wires_dict['yb'][this_wire]

        p0 = np.array([xa, ya])
        p1 = np.array([xb, yb])

        xa, ya = wires_dict['xa'][that_wire], wires_dict['ya'][that_wire]
        xb, yb = wires_dict['xb'][that_wire], wires_dict['yb'][that_wire]

        p2 = np.array([xa, ya])
        p3 = np.array([xb, yb])

        # Find junctions
        J = find_segment_intersection(p0, p1, p2, p3)

        if J is not False:
        # Save coordinates
            xi.append(J[0])
            yi.append(J[1])
            # Save node indices for every edge
            edge_list.append([this_wire, that_wire])

    # Save centres coordinates and edge list to dict
    # if there are junctions
    if len(edge_list) is not 0:
        wires_dict['number_of_junctions'] = len(edge_list)
        wires_dict['xi'] = np.asarray(xi)
        wires_dict['yi'] = np.asarray(yi)
        wires_dict['edge_list'] = np.asarray(edge_list)
        logging.info('Finished detecting junctions')
        return wires_dict
    
    raise Exception('There are no junctions in this network')


def generate_adj_matrix(wires_dict):
    """
generate_wires_distribution    This function will produce adjaceny matrix of 
    the physical network

    Parameters
    ----------

    wires_dict: dict
        a dictionary with all the wires position and junctions/intersection 
        positions.

    Returns
    ------- 
    wires_dict: dict
        The same dictionary with added key:value pairs adjacency matrix 
    """

    # Create array -- maybe use sparse matrix?
    adj_matrix_shape = (wires_dict['number_of_wires'], wires_dict['number_of_wires'])
    adj_matrix = np.zeros(adj_matrix_shape, dtype=np.float32)
    adj_matrix[wires_dict['edge_list'].astype(np.int32)[:, 0], wires_dict['edge_list'].astype(np.int32)[:, 1]] = 1.0
    
    # Make the matrix symmetric
    adj_matrix = adj_matrix + adj_matrix.T

    wires_dict['adj_matrix'] = adj_matrix

    return wires_dict



def generate_graph(wires_dict):
    """
    This function will produce a networkx graph.

    Parameters
    ----------

    wires_dict: dict
        a dictionary with all the wires position and junctions/intersection 
        positions.

    Returns
    ------- 
    wires_dict: dict
        The same dictionary with added key:value pairs networkx graph object.
    """


    # Create graph - this is going to be a memory pig for large matrices
    wires_dict = generate_adj_matrix(wires_dict)
    G = nx.from_numpy_matrix(np.matrix(wires_dict['adj_matrix']))

    wires_dict['G'] = G

    return wires_dict

def analyse_network(wires_dict):
    """
    This function will calculate network statistics for the graph

    Parameters
    ----------

    wires_dict: dict
        a dictionary with all the wires position and junctions/intersection 
        positions.

    Returns
    ------- 
    wires_dict: dict
        The same dictionary with added key:value pairs including network statistics
    """
    
    graph = wires_dict['G']
    
    #diameter of graph
    wires_dict['diameter'] = nx.diameter(graph)
    
    #Shortest path length
    wires_dict['shortpath'] = min(nx.shortest_path_length(graph))
    
    #characteristics path length
    wires_dict['charpath'] = nx.average_shortest_path_length(graph)
    
    #density
    wires_dict['density'] = nx.density(graph)
    
    #circuit rank
    wires_dict['circuit_rank']  = nx.number_of_edges(graph) - nx.number_of_nodes(graph) + nx.number_connected_components(graph)
    
    #average node degree
    wires_dict['avg_nd'] = nx.number_of_edges(graph)*2.0/nx.number_of_nodes(graph)
    
    #standard deviation of node degree
    degrees = np.sum(wires_dict['adj_matrix'], axis=0) #sums along columns
    wires_dict['std_nd'] = np.std(degrees)
    
    return wires_dict



def check_connectedness(wires_dict):
    """
    This function will check is the graph is connected.
    If it is not connected:

    (1) add new junctions 
    (2) update centre coordinates and orientation of one of the nanowires
    (3) something else I haven't thought of yet ...

    """
    
    G = wires_dict['G']

    if not nx.is_connected(G):
        nc = nx.number_connected_components(G)

        logging.warning("This graph has %4d connected components", nc) 
        return False

    return True


def remove_key(wires_dict, key='G'):
    """
    This removes a key:value pair of the dict without 
    altering the original dictionary
    """
    temp_dict = dict(wires_dict)
    del temp_dict[key]
    return temp_dict


def reconnect_graph(wires_dict,scale_num_wires,scale_dispersion,scale_length, falloff_ind_wires=1.5, falloff_ind_length=1.5, falloff_ind_disp=1.5, distro='gamma', loopmax = 10):
    '''
        reconnect graph by increasing length by a set factor each itteration
        more work probably needed such that length and density does not get out of hand.
    '''
    distros = {'gamma': 1, 'uniform':2}

    number_of_wires =int(wires_dict['number_of_wires']*scale_num_wires)
    avg_length = int(wires_dict['avg_length']*scale_length)
    dispersion = int(wires_dict['dispersion']*scale_dispersion)
    length_x = wires_dict['length_x']
    length_y = wires_dict['length_x']
    this_seed = wires_dict['this_seed']
    
    for i in range(loopmax):
        if check_connectedness(wires_dict):
            return wires_dict
        
        number_of_wires =int(wires_dict['number_of_wires']*scale_num_wires)
        avg_length = int(wires_dict['avg_length']*scale_length)
        dispersion = int(wires_dict['dispersion']*scale_dispersion)
        if distros[distro] == 1:
            wires_dict = generate_wires_distribution(number_of_wires,avg_length,dispersion,length_x,length_y,this_seed)
        wires_dict = detect_junctions (wires_dict)
        wires_dict = generate_graph(wires_dict)
        scale_num_wires, scale_length, scale_dispersion= scale_num_wires**falloff_ind_wires, scale_length**falloff_ind_length, scale_dispersion**falloff_ind_disp
    
    return wires_dict



def save_obj(obj, name):
    """
    Save dictionary.
    Might fail in Windows.
    Should be able to give abs path. 
    Currently assuming we are calling from ~/connectivity
    """
    with open('connectivity_data/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    """
    Load object
    """
    with open('connectivity_data/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)



def export_to_matlab(wires_dict, filename=None, save_pkl=False, folder = 'connectivity_data'):
    """
    This exports the dictionary into a matlab file.
    File name convention is as follows:
    Example:
        2016-09-08-153543_asn_nw_02048_nj_11469_seed_042_avl_28.00_disp_10.00_gns_5.00_cdisp_10.00
    timestamp: now
    asn  : atomic switch network
    nw   : number of wires (vertices) 
    nj   : number_of_junctions (edges)
    seed : seed used so we can reproduce the distribution of wire centres
    avl  : average nano wire length
    disp : dispersion of nanowires lengths
    gns  : gennorm shape parameter (beta)
    cdisp: centroid dispersion (scale of gennorm_shape) 
    """
    # Generate a meaningful name

    if filename is None:
        timestamp = time.strftime("%Y-%m-%d-%H%M%S")
        nw   = wires_dict['number_of_wires']
        nj   = wires_dict['number_of_junctions']
        avl  = wires_dict['avg_length']
        disp = wires_dict['dispersion']
        seed = wires_dict['this_seed'] 
        Lx = wires_dict['length_x']
        Ly = wires_dict['length_y']
        if wires_dict['oldNameConvention'] is True:
            gns = wires_dict['gennorm_shape']
            cdisp = wires_dict['centroid_dispersion']
            pars_values = '_asn_nw_%05d_nj_%05d_seed_%03d_avl_%05.2f_disp_%05.2f_gns_%05.2f_cdisp_%05.2f' % (nw, nj, seed, avl, disp, gns, cdisp)
            filename = timestamp + pars_values
        else: 
            filename = 'asn_nw_%05d_nj_%05d_seed_%03d_avl_%05.2f_disp_%05.2f_lx_%05.2f_ly_%05.2f' % (nw, nj, seed, avl, disp, Lx, Ly)
        
        

    # Save the dictionary for later use in python
    if save_pkl:
    	save_obj(wires_dict, filename)

    # Remove key:value that savemat does not understand
    temp_dict = remove_key(wires_dict, 'G')

    # Save to mat format
    import os
    
    #create directory to store file in if it doesn't exist
    if not os.path.exists(folder):
        os.makedirs(folder)

    pathfile = os.path.join(folder, filename + '.mat')
    print('Saved to: ', pathfile)
    savemat(pathfile, temp_dict)
    
    return



def draw_wires(ax, wires_dict):
    """
    Draw wires on a given set of axes.
    
    Wires outside the domain are light gray dashed lines. 
    Wires inside the domain are light gray solid lines. 
    The centre of the wires is marked with a red 'o' marker. 
    
    ax -- matplotlib axes to draw needle symbol
    wires_dict  -- dictionary output from generate_distribution
    """    
    from matplotlib.lines import Line2D
    
    # Make local variables
    xa, ya = wires_dict['xa'], wires_dict['ya']
    xb, yb = wires_dict['xb'], wires_dict['yb']
    xc, yc = wires_dict['xc'], wires_dict['yc']

    for this_wire in range(wires_dict['number_of_wires']):
        if wires_dict['outside'][this_wire]:
            line = [Line2D([xa[this_wire],xb[this_wire]],[ya[this_wire],yb[this_wire]], color='k', ls='--', alpha=0.2)] 
                   #Line2D([xc[this_wire]],[yc[this_wire]],  color='k', marker='o', ms=4, alpha=0.1)] 
        else:   
            line = [Line2D([xa[this_wire],xb[this_wire]],[ya[this_wire],yb[this_wire]], color=(0.42, 0.42, 0.42)),
                    Line2D([xc[this_wire]],[yc[this_wire]], color='r', marker='o', ms=2, alpha=0.77)] 
        for l in line: 
            ax.add_line(l)

    return ax


def draw_junctions(ax, wires_dict):
    """
    Draw the circles at the junctions
    """
    from matplotlib.lines import Line2D
    
    xi, yi = wires_dict['xi'], wires_dict['yi']

    for this_junction in range(wires_dict['number_of_junctions']):
        line = [Line2D([xi[this_junction]],[yi[this_junction]], color='b', marker='o', ms=3, alpha=0.77)]
        for l in line: 
            ax.add_line(l)
    return ax


if __name__ == '__main__':
    # dummy flag to avoid plotting network and importing matplotlib & co
    plot_figures = False

    # Define some parameters
    number_of_wires = 100 
    wire_av_length  = 50.0  # micrometres
    wire_dispersion = 50.0  # micrometres
    this_seed = 42

    wires_dict = generate_wires_distribution(number_of_wires = number_of_wires,
                                       wire_av_length = wire_av_length,
                                       wire_dispersion = wire_dispersion,
                                       this_seed = this_seed)

    detect_junctions(wires_dict)
    # Generate and check graph 
    generate_graph(wires_dict)
    if not check_connectedness(wires_dict):
        logging.warning("Will select the largest connected component.")
        wires_dict = select_largest_component(wires_dict)

    logging.info("The graph is connected. Will save it to mat file.")
    export_to_matlab(wires_dict)

    
    if plot_figures:

        from matplotlib.lines import Line2D
        from matplotlib.patches import Rectangle
        import matplotlib.pyplot as plt

        # Plot pretty pictures of what we just did
        fig, ax = plt.subplots()
        fig.set_size_inches(5,5)

        Lx = wires_dict['length_x']
        Ly = wires_dict['length_y']

        ax.add_patch(Rectangle((0,0), Lx, Ly, color=(1.0, 0.918, 0.0), alpha=0.77))     
        ax = draw_wires(ax, wires_dict)
        ax = draw_junctions(ax, wires_dict)
        ax.set_aspect(1) # set aspect ratio to 1
        ax.set_xlabel(r'x [$\mu$m]')
        ax.set_ylabel(r'y [$\mu$m]')
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.axis([-.1*Lx,1.1*Lx,-.1*Lx,1.1*Lx]) # add some space around the unit square
        ax.set_title('Nanowires distribution')
        ax.grid()
        plt.show()



