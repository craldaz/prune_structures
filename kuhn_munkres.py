#!/usr/bin/env python

import hungarian
import numpy as np
from collections import Counter 

def kabsch_rmsd(traj1, traj2):
    # Subtract COM from each frame
    traj1 = np.array(traj1)
    traj2 = np.array(traj2)
    traj1 -= sum(traj1) / len(traj1) 
    traj2 -= sum(traj2) / len(traj2) 
    
    # Compute covariance matrix
    C = np.dot(np.transpose(traj1), traj2)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
       S[-1] = -S[-1]
       V[:, -1] = -V[:, -1]
 
    # Compute rotation matrix
    U = np.dot(V, W)
    traj1 = np.dot(traj1, U)

    # Compute rmsd from https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
    rmsd = 0.0
    for v, w in zip(traj1, traj2):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(len(traj1[0]))])
    rmsd = np.sqrt(rmsd / len(traj1))
    return rmsd

def sort(labels, coords):
    #List of indices to be sorted
    idx = sorted(range(len(labels)),key=lambda x:labels[x]) 
    labels_sort = [labels[i] for i in idx]
    coords_sort = [list(coords[i,:]) for i in idx]
    return labels_sort, coords_sort

def parse_atoms(labels, coords, atom):
   #Returns a set of coordinates corresponding to parsed atom
   atom_coords = [coords[i] for i in range(len(labels)) if labels[i] == atom]
   return atom_coords

def transform_atoms(coords, swap, reflect, idx):
   #Returns coordinates after transforming specific atoms
   new_coords = [x[:] for x in coords]
   new_coords[i][0] = [coords[i][swap[0]]*reflect[0] for i in idx]
   new_coords[i][1] = [coords[i][swap[1]]*reflect[1] for i in idx]
   new_coords[i][2] = [coords[i][swap[2]]*reflect[2] for i in idx]
   return new_coords

def transform_coords(coords, swap, reflect):
   #Returns the transformed coordinates
   new_coords = [[coords[i][swap[0]]*reflect[0],coords[i][swap[1]]*reflect[1],coords[i][swap[2]]*reflect[2]] for i in range(len(coords))]
   return new_coords

def permute_atoms(coords, permutation, idx):
   #Returns the coordinates after permuting just the specified atom
   new_coords = coords[:]
   for i in range(len(permutation)):
      j = idx[permutation[i]]
      k = idx[i]
      new_coords[k] = coords[j] 
   return new_coords

def parse_indices(labels, atom):
   #Return indices of labels and atoms
   indices = [i for i in range(len(labels)) if labels[i] == atom]
   return indices

def build_hashtable(labels,coords, Uniq, x):
    Coords = {}
    Indices = {}
    for i in range(len(Uniq)):
        Coords[Uniq[i]]  = parse_atoms(labels, coords, str(Uniq[i]))
        Indices[Uniq[i]] = parse_indices(labels, str(Uniq[i]))
    return Coords, Indices

def kuhn_munkres(traj1,traj2):
    #Atom labels/coords/# atoms 
    NA_a = traj1.n_atoms
    NA_b = traj2.n_atoms
    a_coords = np.reshape(traj1.xyz, (traj1.n_atoms,3))
    b_coords = np.reshape(traj2.xyz, (traj1.n_atoms,3))
    a_labels = [str(traj1.top.atom(i))[2:] for i in range(traj1.n_atoms)]    
    b_labels = [str(traj2.top.atom(i))[2:] for i in range(traj2.n_atoms)]    
    
    #Compute RMSD Standard way
    InitRMSD_unsorted = kabsch_rmsd(a_coords,b_coords)
    #print('Unsorted atom RMSD: %.4f' % InitRMSD_unsorted)
    
    #Sort the atom labels and coords
    a_labels, a_coords = sort(a_labels,a_coords)
    b_labels, b_coords = sort(b_labels,b_coords)

    #Count number of unique atoms
    Uniq = Counter(a_labels).keys()
    n_types = len(Counter(a_labels).values())
  
    Perm = {}
    for i in range(len(Uniq)):
        Perm[Uniq[i]] = 'perm_' + Uniq[i]
        Perm[Uniq[i]] = []
     
    #Generate hastables for traj1 and traj2 and remove COM
    a_Coords, a_Indices = build_hashtable(a_labels, a_coords, Uniq, 'a_')
    b_Coords, b_Indices = build_hashtable(b_labels, b_coords, Uniq, 'b_')
    A = np.array(a_Coords[Uniq[0]])
    B = np.array(b_Coords[Uniq[0]])
    A -= sum(A)/len(A)
    B -= sum(B)/len(B)
     
    #Define swap and reflections for initial B (coordinates of traj2) and apply
    swaps = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    reflects = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1),(-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1)]
    B_transformed = [[transform_coords(B,i,j),i,j] for i in swaps for j in reflects]
    
    #Performs the munkres algorithm on each set of transformed coordinates
    rmsds = []
    for i in range(len(B_transformed)):
        l = 0
        cost_matrix = np.array([[np.linalg.norm(a - b) for b in B_transformed[i][0]] for a in A])
        LAP = hungarian.lap(cost_matrix)
        
        Perm[Uniq[l]] = []
        for j in range(len(LAP[0])):
           Perm[Uniq[l]] += [(j,LAP[0][j])]
        Perm[Uniq[l]] = sorted( Perm[Uniq[l]], key = lambda x: x[0])
        Perm[Uniq[l]] = [x[1] for x in Perm[Uniq[l]]]
        b_perm = permute_atoms(b_coords, Perm[Uniq[l]], b_Indices[Uniq[l]])
        b_trans = transform_coords(b_perm, B_transformed[i][1], B_transformed[i][2])
        while l < n_types:
            if l > 0:
               b_Coords[Uniq[l]] = parse_atoms(b_labels, b_final, Uniq[l])
            else:
               b_Coords[Uniq[l]] = parse_atoms(b_labels, b_trans, Uniq[l])
            cost_matrix = np.array([[np.linalg.norm(a- b) for b in np.array(b_Coords[Uniq[l]])] for a in np.array(a_Coords[Uniq[l]])])
            LAP = hungarian.lap(cost_matrix)
            Perm[Uniq[l]] = []
            for k in range(len(LAP[0])):
               Perm[Uniq[l]] += [(k,LAP[0][k])]
            Perm[Uniq[l]] = sorted( Perm[Uniq[l]], key = lambda x: x[0])
            Perm[Uniq[l]] = [x[1] for x in Perm[Uniq[l]]]
            b_final = permute_atoms(b_trans, Perm[Uniq[l]], b_Indices[Uniq[l]])
            b_trans = b_final
            l += 1
            q = l - 1 
            rmsds.append([kabsch_rmsd(a_coords, b_final), B_transformed[i][1], B_transformed[i][2], b_final])
            rmsds = sorted(rmsds, key = lambda x: x[0])
    #print('Hungarian Algorithm RMSD: %.4f' % rmsds[0][0])
    return rmsds[0][0]
