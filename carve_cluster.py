#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 15:56:44 2018
Modified on Tues Dec 04 11:11:13 2018

@author: ernesto
"""
##############################################################################
# This code reads a user-INPUT xyz file with a central metallic center
# and extracts clusters from the center out based on the atom name (solvent)
# for example O for OH or H2O.

# Modified to include other species in the solvated cluster, i.e. al-n(h20)
# with na at different distances.
##############################################################################
import sys, os 
import numpy as np
import pandas as pd
import math

#-----------------------------------------------------------------------------#
#Some directives to get data from xyz file

print('In this code you need to INPUT: ' \
      'Name of the xyz file as sys.argv[1] ' \
      'Number of solvent molecules around the metal center'\
      'Label of the central atom')
#print('Input the name of the xyz file to be carved out \n')
infile = open(sys.argv[1], 'r')
#print('How many molecules you want the cluster to be formed of \n')
mol_number = sys.argv[2]
l_atom = sys.argv[3]        #Label of the central atom
extra_atom = sys.argv[4]    #Add an extra atom label to the xyz / yes or no
e_atom = sys.argv[5]        #Label of the extra atom to add in each cluster or type no

print(f'Working on file {infile.name}')

#-----------------------------------------------------------------------------#
# Initializing input file reading and outputting file
infile_lines = infile.readlines()
if infile_lines:
    if infile_lines[0]:
        N = int(infile_lines[0].split()[0])
else:
    print(f"The file {infile.name} doesn't have content, or there is a " \
           "problem with the initial N in file.xyz")

snaps = int(len(infile_lines)/(N+2))
print(f'There are {snaps} snapshots in the {infile.name} file')

outf_name = 'cluster.xyz'
try:
    os.remove(outf_name)
    print(f'Previous file {outf_name} was removed')
except OSError:
    print(f"Output file didn't exist before")
    pass
os.system(f'mkdir cluster_N_{mol_number}')
    
outfile_b = open(outf_name, 'ab')

#-----------------------------------------------------------------------------#
# Loop over snapshots, over element types to form molecules based on number of
# molecules requested.

for i in range(snaps):
    snap_xyz = []
    d = []
    for j in range(N+2):
        if j == 0 or j == 1:
            continue
        snap_xyz.append(infile_lines[i*(N+2)+j].split())                        
    dataf_xyz = pd.DataFrame(data = snap_xyz, columns = ['Label', 'X', 'Y', 'Z'], dtype = float)
    print(dataf_xyz[dataf_xyz['Label'] == 'Li'].values)
    l = np.zeros((len(dataf_xyz[dataf_xyz['Label']=='O']), 2), dtype=float)
    for j in range(len(dataf_xyz[dataf_xyz['Label']=='O'])):
        distance = math.sqrt(math.pow(dataf_xyz[dataf_xyz['Label'] == 'O'].iloc[j]['X']-dataf_xyz[dataf_xyz['Label'] == l_atom]['X'], 2) + \
                   math.pow(dataf_xyz[dataf_xyz['Label'] == 'O'].iloc[j]['Y']-dataf_xyz[dataf_xyz['Label'] == l_atom]['Y'], 2) + \
                   math.pow(dataf_xyz[dataf_xyz['Label'] == 'O'].iloc[j]['Z']-dataf_xyz[dataf_xyz['Label'] == l_atom]['Z'], 2)  )
        l[j][0] = dataf_xyz[dataf_xyz['Label']=='O'].index[j]
        l[j][1] = distance
#   print(np.shape(l))
    d_array = l[l[:,1].argsort()]
    p = np.zeros((len(dataf_xyz[dataf_xyz['Label']=='H']), 2), dtype=float)
    cluster = pd.DataFrame(columns = ['Label','X', 'Y', 'Z'])
    cluster = cluster.append(dataf_xyz[dataf_xyz['Label'] == l_atom])
    for k in range(int(mol_number)):
        for m in range(len(dataf_xyz[dataf_xyz['Label']=='H'])):
            o_h_distance = math.sqrt(math.pow(dataf_xyz[dataf_xyz['Label'] == 'H'].iloc[m]['X']-dataf_xyz.iloc[int(d_array[:][k][0])]['X'], 2) + \
                           math.pow(dataf_xyz[dataf_xyz['Label'] == 'H'].iloc[m]['Y']-dataf_xyz.iloc[int(d_array[:][k][0])]['Y'], 2) + \
                           math.pow(dataf_xyz[dataf_xyz['Label'] == 'H'].iloc[m]['Z']-dataf_xyz.iloc[int(d_array[:][k][0])]['Z'], 2)  )
            p[m][0] = dataf_xyz[dataf_xyz['Label']=='H'].index[m]
            p[m][1] = o_h_distance
#        print(np.shape(p))
        d_array_H = p[p[:,1].argsort()]
        s = []
        for o in range(len(d_array_H[d_array_H[:,1] < 1.2][:,0])):
            s.append(int(d_array_H[d_array_H[:,1] < 1.2][:,0][o]))
        b = int(d_array[:][k][0])
#        print(type(dataf_xyz.iloc[[b]]), '\n', type(dataf_xyz.iloc[s]))
#        print(dataf_xyz.iloc[[b]], '\n', dataf_xyz.iloc[s])
        cluster = cluster.append(dataf_xyz.iloc[[b]], ignore_index = True)
        cluster = cluster.append(dataf_xyz.iloc[s], ignore_index = True)
# -----------------------------------------------------------------------------#
# Imprime data the l_atom, O, and H to form requested number of molecules in the cluster
#    print(type(cluster))
# Adding an extra atom label with its coordinates
    if extra_atom == 'yes':
        N_atom = len(cluster) + 1
        with open(outf_name, 'a+') as outfile_s:
            outfile_s.write(f'{N_atom} \n\n')
        with open(outf_name, 'ab') as outfile_b:
            np.savetxt(outfile_b, cluster.values, fmt='%s')
            np.savetxt(outfile_b, dataf_xyz[dataf_xyz['Label'] == e_atom].values, fmt='%s')
    else:
        N_atom = len(cluster)
        with open(outf_name, 'a+') as outfile_s:
            outfile_s.write(f'{N_atom} \n\n')
        with open(outf_name, 'ab') as outfile_b:
            np.savetxt(outfile_b, cluster.values, fmt='%s')

infile.close()
os.system(f'mv {outf_name} cluster_N_{mol_number}/.')