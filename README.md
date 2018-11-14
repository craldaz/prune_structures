# prune_structures

print "Usage: python prune.py -e [energy_cutoff] -n [max_clusters -r [rmsd_threshold] -a [aligned_atoms] --smi [reference smi] -m [cluster_mode]"
          -e   Energy cutoff for conformers. Positive flaot (kCal/mol)   # look for second line 
          -n   Max number of clusters to be kept. Positive integer.      #
          -r   RMSD threshold for conformers within a cluster
          -a   List of aligned atoms IDs (0-indexed) seperated by comma
               eg: "0, 3, 5"
          --smi reference SMILE string to prune wrong molecules
          -m   Clustering mode
               0: by RMSD (Hungarian Algorithm) #This is very slow and should be used with under 10-15 structures
               1: by dist spectrum

