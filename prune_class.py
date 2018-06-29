from generate import readXyzEnergy
from generate import printLine
import subprocess
import sys
import os
import getopt
import glob
from time import sleep
import pybel
import mdtraj as md
import numpy as np
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform
from scipy.spatial import distance
import kuhn_munkres as km

#Some constants
kCaltoAu = 0.00159360109742136
#Each entry of the distance matrix is based on the rmsd of two molecules
def clustering(traj, maxClusters, thresh_R):
  distances=np.empty((traj.n_frames,traj.n_frames))
  for i in range(traj.n_frames):
      for j in range(traj.n_frames):
          if md.rmsd(traj[i],traj[j],0) <= 0.001:
              distances[i,j] = 0.0
          else:
              distances[i,j] = km.kuhn_munkres(traj[i], traj[j])
  mean_rmsd = np.mean(distances[np.nonzero(distances)])
  min_rmsd = np.min(distances[np.nonzero(distances)])
  max_rmsd = np.max(distances[np.nonzero(distances)])
  print('Min pairwise rmsd: %f nm' % np.min(min_rmsd))
  print('Max pairwise rmsd: %f nm' % np.max(max_rmsd))
  print('Mean pairwise rmsd: %f nm' % np.max(mean_rmsd))
  reduced_distances = squareform(distances, checks=False)
  linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='average')
  flat_cluster=scipy.cluster.hierarchy.fcluster(linkage,thresh_R,criterion='distance')
  return flat_cluster

def checkAndRemoveFile(nameWildcard):
    print 'Checking existence of old file(s)',nameWildcard
    fileList=glob.glob(nameWildcard)
    if len(fileList) > 0:
       print 'Old file(s) ',nameWildcard,'exist in currect directory'
       for f in fileList:
          subprocess.call(['rm', f])
       print 'Removed old file(s) to avoid confusion'

class prune:
    def __init__(self):
        self.workDir_ = ""
        self.eThresh_ = 20
        self.maxClusters_ = 50
        self.rmsdThresh_ = 0.2
        self.alignedAtoms_ = []
        self.nAtoms_ = 0
        self.numSortedMin_ = 0
        self.numEffectMin_ = 0
        self.mdXyzFileName_ = 'all.xyz'
        self.nameEnergyList_ = []
        self.clusterMode_ = 0
        self.t_aligned_ = None
        self.step_ = 0
        self.refSmi_ = ""
        self.spectrumDifference_ = None

    def printHeader(self):
        print """ 
        ##########################
        # Starting AUTO NEB run  #
        ##########################
        Prune the conformers by energy and clustering... """       

    def usage(self):
        print "Usage: python prune.py -e [energy_cutoff] -n [max_clusters -r [rmsd_threshold] -a [aligned_atoms] --smi [reference smi] -m [cluster_mode]"
        print ''' 
                  -e   Energy cutoff for conformers. Positive flaot (kCal/mol)
                  -n   Max number of clusters to be kept. Positive integer.
                  -r   RMSD threshold for conformers within a cluster
                  -a   List of aligned atoms IDs (0-indexed) seperated by comma
                       eg: "0, 3, 5"
                  --smi reference SMILE string to prune wrong molecules
                  -m   Clustering mode
                       0: by RMSD (hungarian)
                       1: by dist spectrum
                       '''

    def readOpt(self):
        print 'Argument List:', str(sys.argv)
        try:
          optlist,args = getopt.getopt(sys.argv[1:],'e:n:r:a:m:',['smi='])
        except getopt.GetoptError as err:
          print str(err)
          self.usage()
          sys.exit(2)   
        for o,a in optlist:
             if o == "-e":
               self.eThresh_ = float(a) 
             elif o == "-n":
               self.maxClusters_ = int(a) 
             elif o == "-r":
               self.rmsdThresh_ = float(a)
             elif o == "-a":
               self.alignedAtoms_ = [int(k) for k in a.split(",")] 
             elif o == "--smi":
               self.refSmi_ = a
             elif o == "-m":
               self.clusterMode_ =int(a)
             else:
               assert False, "unhandled option"

    def printSanityCheck(self,term,val,append=""):
        print term.ljust(40),val,append

    def sanityCheck(self):
        self.printHeader()
        self.workDir_ = subprocess.Popen(['pwd'],
                         stdout=subprocess.PIPE).communicate()[0].strip('\n')
         
        if self.eThresh_ < 0:
             assert False, "Energy threshold must be a positive float number"
        if self.maxClusters_ <0:
             assert False, "Maximum number of clusters must be a positive integer"
        if self.rmsdThresh_ < 0:
             assert False, "Threshold of RMSD within a cluster must a positive float"
        try:
            f=open('conf0002.xyz')
        except IOError as e:
            print 'Conformer structure file conf0.xyz does not exist in current directory'
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
        self.nAtoms_ = len(f.readlines())-2
            
        if self.alignedAtoms_ == []:
             print '''WARNING: Atoms to be aligned is not specified. 
                      By default, all atoms will be aligned to reference
                   '''
             self.alignedAtoms_ = [ int(k) for k in range(0,self.nAtoms_) ]
        if self.refSmi_ !="":
             self.printSanityCheck('Molecules will be checked against SMILE string:',self.refSmi_)
        else:
             print 'WARNING: No reference SMILE string provided. Molecules will not be checked.'
        
        if self.clusterMode_ == 0:
             self.printSanityCheck('Molecules will be clustered based on RMSD Computed byt Hungarian Algorithm',"")
        elif self.clusterMode_ ==1:
             self.printSanityCheck('Molecules will be clustered based on distance spectrum',"")
        else:
             print "ERROR: unhandled type of clustering mode"
             sys.exit(2)

        self.printSanityCheck('My work directory is ',self.workDir_)
        self.printSanityCheck('Energy cutoff for conformers:',self.eThresh_,'kCal/mol')
        self.printSanityCheck('Max number of clusters:',self.maxClusters_)
        self.printSanityCheck('RMSD threshold for each cluster:', self.rmsdThresh_)
        self.printSanityCheck('Following atoms will be aligned:',self.alignedAtoms_)
        printLine()
        self.step_+=1
        

####################
#  Read Conformers #
####################
    def readConformers(self):
        print self.step_,'. Reading all energetically sorted conformers...'
        cmd = 'ls conf*.xyz'
        proc1 = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        cmd = ['wc','-l']
        proc2 = subprocess.Popen(cmd,stdin=proc1.stdout,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        self.numSortedMin_= int(proc2.communicate()[0])
        if self.numSortedMin_ < 1:
           print 'No minimized conformer found in current directory. Nothing to do.'
           print 'Please check the unconverged minimization runs and proceed again.'
           return
        
        print 'Number of sorted, minimized conformers:',self.numSortedMin_
        
        self.nameEnergyList_=[None]*self.numSortedMin_
        for i in range(2,self.numSortedMin_):
           #name='conf'+str(i)+'.xyz'
           name='conf'+'{:04}'.format(i)+'.xyz'
           energy=readXyzEnergy(name)
           self.nameEnergyList_[i]=[name,energy]
        printLine()
        self.step_ +=1

####################
#  Prune by energy #
####################
    def pruneByEnergy(self):
        print self.step_,'. Prune conformers by energy'
        print 'Conformers with energy higher than ', self.eThresh_, 'kCal/mol are removed'
        print 'Remaining conformers are saved to all.xyz'
        refE = self.nameEnergyList_[0][1]
        self.numEffectMin_ = 0
        for i in range(0,len(self.nameEnergyList_)):
           deltaE = (self.nameEnergyList_[i][1] - refE)/kCaltoAu
           self.nameEnergyList_[i][1] = deltaE
           if deltaE > self.eThresh_:
              break
           else:
              self.numEffectMin_ +=1
        print 'After pruning by energy, number of effective conformers:', self.numEffectMin_
        del self.nameEnergyList_[self.numEffectMin_:]

        printLine()
        self.step_ +=1

##################################################
# Remove wrong molecules by checking fingerprint #
##################################################
    def fingerPrintChecking(self):
       if self.refSmi_ !="":
           print self.step_,'. Delete structures not matching the SMILES string'
           refMol=pybel.readstring("smi",self.refSmi_)
           refFp =refMol.calcfp()
           mols = [pybel.readfile("xyz",x[0]).next() for x in self.nameEnergyList_]      
           fps = [x.calcfp() for x in mols]
           TanimotoCoeff = [ refFp | fps[i] for i in range(0,len(fps)) ]
           removeIndices = [ i for i in range(0,len(fps) ) if TanimotoCoeff[i] < 1 ]
           removeNames = [ self.nameEnergyList_[i][0] for i in removeIndices ]
           print "The following xyz files are removed because the molecule structure doesn't match the reference SMILE string:"
           print ("\n".join(removeNames))
           # Becareful we should remove the indices in descending order inorder not to change the labeling
           for i in reversed(removeIndices):
                 del self.nameEnergyList_[i]
           self.numEffectMin_ = len(self.nameEnergyList_)
           print "Updated number of effective conformers:", self.numEffectMin_
           printLine()
           self.step_ +=1
##############################################
#  Remove duplicated conformers by checking  #
#  eigen values of geometric distance matrix #
##############################################
    '''Theoretical background
       Geometric distance matrix:
       is the distance between each pair of atoms in this molecule. 
       The eigenvalues of the distance matrix of a conformer is called the distance spectrum.
       Difference conformers should have unique distance spectrum
       Because of the symmetry or different labeling of equivalent atoms, we could have
       duplicated conformers whose distance spectrum are the same. 
       Based on these we can remove duplicatted conformers.
    '''
    def geometricSpectrumChecking(self):
        print self.step_,'. Prune conformers by geometry spectrum'
        self.numEffectMin_ =  len(self.nameEnergyList_)
        distSpectrums = np.empty([self.numEffectMin_,self.nAtoms_])
        for i in range(0,self.numEffectMin_):
            name = self.nameEnergyList_[i][0]
            mol = pybel.readfile('xyz',name).next()
            coords = [iatom.coords for iatom in mol.atoms]
            distMatrix= squareform(distance.pdist(np.asarray(coords).reshape(-1,3)))
            [distSpect,egv] = np.linalg.eigh(distMatrix)
            distSpectrums[i][:]=distSpect[:]
        self.spectrumDifference_=squareform(distance.pdist(distSpectrums))
        diffThresh = 0.01
        removeIndices = []
        #when two conformers are duplicated, keep the one with smaller index
        for i in range(self.numEffectMin_-2, -1, -1):
            for j in range(self.numEffectMin_-1, i, -1):
               diff = np.sqrt(self.spectrumDifference_[i][j]/self.nAtoms_)
               if diff < diffThresh:
                   print self.nameEnergyList_[i][0],' and ',self.nameEnergyList_[j][0], ' are likely to be duplicated conformers. Dist Spectrum diff=',diff
                   if not (j in removeIndices):
                      removeIndices.append(j)
        removeIndices = sorted(removeIndices, reverse=True)
        for i in removeIndices:
            name = self.nameEnergyList_[i][0]
            del self.nameEnergyList_[i]
            print 'Remove conformer ', name

        self.numEffectMin_ =  len(self.nameEnergyList_)
        print "Updated number of effective conformers:", self.numEffectMin_
        printLine()
        self.step_ +=1
#####################
#  Align structures #
#####################
    def align(self):
        print self.step_, 'Align all effective conformers to the lowest-energy conformer'
        checkAndRemoveFile(self.mdXyzFileName_)
        print 'Wring out effective conformers to file:',self.mdXyzFileName_
        mdXyzFile=open(self.mdXyzFileName_,'a')
        for i in range(0,self.numEffectMin_):
           name=self.nameEnergyList_[i][0]
           subprocess.call(['cat',name],stdout=mdXyzFile)
        mdXyzFile.close()
        #Get Topology file to be used in mdtraj
        topXyzName = self.nameEnergyList_[0][0]
        topPdbName = topXyzName.replace('xyz','pdb')
        checkAndRemoveFile(topPdbName)
        mol = pybel.readfile('xyz',topXyzName).next()
        output = pybel.Outputfile('pdb',topPdbName, overwrite=True)
        output.write(mol)
        output.close()
        
        #Load conformers into mdtraj and analyze
        t=md.load(self.mdXyzFileName_,top=topPdbName)
        self.t_aligned_ =t.superpose(t,frame=0,atom_indices=self.alignedAtoms_)
        checkAndRemoveFile('all_aligned.xyz')
        self.t_aligned_.save_xyz('all_aligned.xyz')
        printLine()
        self.step_ +=1
#######################
# Prune by clustering #
#######################
#Each entry of distances matrix is based on comparison of dist spectrums of two molecules
    def clustering2(self, traj, maxClusters, thresh_R):
        distances=np.empty((traj.n_frames,traj.n_frames))
        for i in range(traj.n_frames):
           for j in range(traj.n_frames):
              distances[i][j] = self.spectrumDifference_[i][j]/self.nAtoms_
        mean_rmsd = np.mean(distances[np.nonzero(distances)])
        min_rmsd = np.min(distances[np.nonzero(distances)])
        max_rmsd = np.max(distances[np.nonzero(distances)])
        print('Min pairwise rmsd: %f nm' % np.min(min_rmsd))
        print('Max pairwise rmsd: %f nm' % np.max(max_rmsd))
        print('Mean pairwise rmsd: %f nm' % np.max(mean_rmsd))
        reduced_distances = squareform(distances, checks=False)
        linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method='average')
        flat_cluster=scipy.cluster.hierarchy.fcluster(linkage,thresh_R,criterion='distance')
        return flat_cluster

    def pruneByClustering(self):
        print self.step_, '. Prune conformers by clustering'
        #Clustering the structures
        if self.clusterMode_ == 0:
          clusterInd=clustering(self.t_aligned_, self.maxClusters_, self.rmsdThresh_)
        else:
          clusterInd=self.clustering2(self.t_aligned_, self.maxClusters_, self.rmsdThresh_)
      
        #Based on the clustering, print the conformers of each cluster
        '''Build a map:
              cluster ID <==> List of frame numbers
        '''
        clusterToFrameDict={}
        for i in range(0,len(clusterInd)):
           myClusterID = clusterInd[i]
           myFrameID  = i
           if myClusterID in clusterToFrameDict:
             clusterToFrameDict[myClusterID].append(i)
           else:
             clusterToFrameDict[myClusterID]=[i]
        print len(clusterToFrameDict), ' clusters are generated'
        print clusterToFrameDict
      
        '''Save the clustered structures to file. 
           1. For each cluster, save the conformer with lowest energy 
              (i.e. the one with smallest frame ID, because originally the 
                conformers are ordered by energy)
           2. Define the energy of each cluster as the energy of the lowest-
              energy former within that cluster
           3. Sort clusters by their energy
           4. Print clusters to file
           
        '''
        clusterEnergyList=[]
        for key in clusterToFrameDict:
            myLowestEnergyConfID = clusterToFrameDict[key][0]
            myLowestEnergy = self.nameEnergyList_[myLowestEnergyConfID][1]
            clusterEnergyList.append([myLowestEnergyConfID,myLowestEnergy])
        clusterEnergyList.sort(key=lambda tup:tup[1])
        digitsToPad=int(np.ceil(np.log10(len(clusterEnergyList))))
        checkAndRemoveFile('cluster*.xyz')
        alignedmols=list(pybel.readfile('xyz','all_aligned.xyz'))
        for i in range(0,len(clusterEnergyList)):
            outputname = 'cluster'+str(i).zfill(digitsToPad)+'.xyz'
            output = open(outputname,'w')
            myFrameID = clusterEnergyList[i][0] 
            #self.t_aligned_[myFrameID].save_xyz(outputname)
            myOBMol=alignedmols[myFrameID].OBMol
            info=str(clusterEnergyList[i][1])+' kcal/mol'
            myOBMol.SetTitle(info)
            pybel.Molecule(myOBMol).write('xyz',outputname,overwrite=True)
            print 'Cluster',i,'conformer id',myFrameID, 'Energy',clusterEnergyList[i][1]

#########################
#      Main Driver      #
#########################
    def run(self):
        self.readOpt()
        self.sanityCheck()
        self.readConformers()
        #self.fingerPrintChecking()            #CALDAZ TURNED OFFF
        self.geometricSpectrumChecking()
        self.pruneByEnergy()
        self.align()
        self.pruneByClustering()
         
if __name__ == '__main__':
    p = prune()
    p.run()
