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
import networkx as NX
import pygraphviz as PG

#Some constants
kCaltoAu = 0.00159360109742136

def demoBarrierMatrix(n):
    b=np.random.random_sample((n,n))*50
    b_symm=(b+b.T)*0.5
    for i in range(0,n):
      b_symm[i][i]=0
    return b_symm
    

def getDistMatrix(traj):
  distances=np.empty((traj.n_frames,traj.n_frames))
  for i in range(traj.n_frames):
     distances[i] = md.rmsd(traj, traj, i)
  return distances

def checkAndRemoveFile(nameWildcard):
    print 'Checking existence of old file(s)',nameWildcard
    fileList=glob.glob(nameWildcard)
    if len(fileList) > 0:
       print 'Old file(s) ',nameWildcard,'exist in currect directory'
       for f in fileList:
          subprocess.call(['rm', f])
       print 'Removed old file(s) to avoid confusion'

def minDist(distMatrix):
    if len(distMatrix) < 2:
       assert False, 'DistMatrix should have dimension >= 2'
    dist= distMatrix[0][1]
    [n,m]= np.shape(distMatrix)
    for i in range(0,n):
       for j in range(i+1,n):
          if distMatrix[i][j] < dist:
              dist = distMatrix[i][j] 
    return dist

def maxDist(distMatrix):
    if len(distMatrix) < 2:
       assert False, 'DistMatrix should have dimension >= 2'
    dist= distMatrix[0][1]
    [n,m]= np.shape(distMatrix)
    for i in range(0,n):
       for j in range(i+1,n):
          if distMatrix[i][j] > dist:
              dist = distMatrix[i][j] 
    return dist

def plotConfNet(nameList, energyList, DistMatrix,BarrierMatrix=[], LengthScale=2.0):
    G = PG.AGraph()
    G.add_nodes_from(range(len(nameList)))
    G.node_attr.update(color='red', style='rounded', shape='box', overlap= 'prism', splines='true')
    MaxEdgeLength = 20
    DistinguishableLength = 0.5
    lenFactor= min(LengthScale/minDist(DistMatrix),MaxEdgeLength/maxDist(DistMatrix)) 
    print 'Edge Length Scaling:',lenFactor
    if BarrierMatrix==[]:
        for i in range(0,len(nameList)):
            for j in range(i+1,len(nameList)):
               length=str(lenFactor*DistMatrix[i][j])
               G.add_edge(i, j, len=length, color = 'blue')
               #G.add_edge(i, j, color = 'blue')
               if lenFactor*DistMatrix[i][j] < DistinguishableLength:
                  print 'WARNING: These conformers seems to be identical and not distinguishable on the graph:', i,j, 'Dist', DistMatrix[i][j]
    else:
        minPenwidth = 0.5
        maxBarrier =max(BarrierMatrix.flatten())
        for i in range(0,len(nameList)):
            for j in range(i+1,len(nameList)):
               length=str(lenFactor*DistMatrix[i][j])
               G.add_edge(i, j, len=length, color = 'blue', penwidth = str(minPenwidth*maxBarrier/BarrierMatrix[i][j]))
               #G.add_edge(i, j, color = 'blue', penwidth = str(minPenwidth*maxBarrier/BarrierMatrix[i][j]))
               if lenFactor*DistMatrix[i][j] < DistinguishableLength:
                  print 'WARNING: These conformers seems to be identical and not distinguishable on the graph:', i,j, 'Dist', DistMatrix[i][j]

    for i in range(0,len(nameList)):
        n=G.get_node(i)
        n.attr['label'] = nameList[i].strip('.xyz')+'\n( '+str(round(energyList[i],2))+' kCal/mol)'
    G.draw('conformers.png', format='png', prog='neato')

class plotConformers:
    def __init__(self):
        self.workDir_ = ""
        self.eThresh_ = 10
        self.alignedAtoms_ = []
        self.nAtoms_ = 0
        self.numSortedMin_ = 0
        self.numEffectMin_ = 0
        self.mdXyzFileName_ = 'all.xyz'
        self.nameEnergyList_ = []
        self.prefix_ = 'conf'
        self.t_aligned_ = None
        self.unitConversion_ = False
        self.demo_ = False
        self.step_ = 0

    def printHeader(self):
        print """ 
        ##########################
        # Starting AUTO NEB run  #
        ##########################
        Prune the conformers by energy and clustering... """       

    def usage(self):
        print "Usage: python plot.py -e [energy_cutoff] -a [aligned_atoms] -p [conformer_name_prefix] -u [energy_unit_conversion]"
        print ''' 
                  -e   Energy cutoff for conformers. Positive flaot (kCal/mol)
                  -a   List of aligned atoms IDs (0-indexed) seperated by comma
                       eg: "0, 3, 5"
                  -p   conformer name prefix. conformers should have names like [prefix][0,1,2,...].xyz
                       default prefix is 'conf'
                  -u   energy unit conversion from Hartree to Kcal/mol.
                  -d   demo mode with barrier matrix
                       '''

    def readOpt(self):
        print 'Argument List:', str(sys.argv)
        try:
          optlist,args = getopt.getopt(sys.argv[1:],'e:a:p:ud')
        except getopt.GetoptError as err:
          print str(err)
          self.usage()
          sys.exit(2)   
        for o,a in optlist:
             if o == "-e":
               self.eThresh_ = float(a) 
             elif o == "-a":
               self.alignedAtoms_ = [int(k) for k in a.split(",")] 
             elif o == '-p':
               self.prefix_=a
             elif o == '-u':
               self.unitConversion_ = True
             elif o == '-d':
               self.demo_ = True
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
        try:
            f=open(self.prefix_+'0.xyz')
        except IOError as e:
            print 'Conformer structure file ',self.prefix_,'0.xyz does not exist in current directory'
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
        
        self.nAtoms_ = len(f.readlines())-2
            
        if self.alignedAtoms_ == []:
             print '''WARNING: Atoms to be aligned is not specified. 
                      By default, all atoms will be aligned to reference
                   '''
             self.alignedAtoms_ = [ int(k) for k in range(0,self.nAtoms_) ]

        self.printSanityCheck('My work directory is ',self.workDir_)
        self.printSanityCheck('Energy cutoff for conformers:',self.eThresh_,'kCal/mol')
        self.printSanityCheck('Following atoms will be aligned:',self.alignedAtoms_)
        printLine()
        self.step_+=1

####################
#  Read Conformers #
####################
    def readConformers(self):
        print self.step_,'. Reading all energetically sorted conformers...'
        cmd = 'ls '+self.prefix_+'*.xyz'
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
        for i in range(0,self.numSortedMin_):
           name=self.prefix_+str(i)+'.xyz'
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
        for i in range(0,self.numSortedMin_):
           if(self.unitConversion_ ==True):
             deltaE = (self.nameEnergyList_[i][1] - refE)/kCaltoAu
             self.nameEnergyList_[i][1] = deltaE
           else:
             deltaE = (self.nameEnergyList_[i][1] - refE)
             self.nameEnergyList_[i][1] = deltaE
             
           if deltaE > self.eThresh_:
              break
           else:
              self.numEffectMin_ = i+1
        print 'After pruning by energy, number of effective conformers:', self.numEffectMin_
        del self.nameEnergyList_[self.numEffectMin_:]

        checkAndRemoveFile(self.mdXyzFileName_)
        print 'Wring out energy pruned conformers to file:',self.mdXyzFileName_
        mdXyzFile=open(self.mdXyzFileName_,'a')
        for i in range(0,self.numEffectMin_):
           name=self.prefix_+str(i)+'.xyz'
           subprocess.call(['cat',name],stdout=mdXyzFile)
        mdXyzFile.close()
        printLine()
        self.step_ +=1

#####################
#  Align structures #
#####################
    def align(self):
        print self.step_, 'Align all conformers to the lowest-energy conformer'
        #Get Topology file to be used in mdtraj
        topPdbName = self.prefix_+'0.pdb'
        checkAndRemoveFile(topPdbName)
        topXyzName = self.prefix_+'0.xyz'
        mol = pybel.readfile('xyz',topXyzName).next()
        output = pybel.Outputfile('pdb',topPdbName, overwrite=True)
        output.write(mol)
        output.close()
        
        #Load conformers into mdtraj and analyze
        t=md.load(self.mdXyzFileName_,top=topPdbName)
        self.t_aligned_ =t.superpose(t,frame=0,atom_indices=self.alignedAtoms_)
        printLine()
        self.step_ +=1
##########################
# Draw conformer network #
##########################
    def plot(self):
        print self.step_, '. Plot conformers as a network'
        DistMatrix=getDistMatrix(self.t_aligned_)
        nameList=[ row[0] for row in self.nameEnergyList_ ]
        energyList=[ row[1] for row in self.nameEnergyList_ ]
        if self.demo_ ==True:
          plotConfNet(nameList,energyList,DistMatrix,demoBarrierMatrix(len(energyList)) )
        else:
          plotConfNet(nameList,energyList,DistMatrix)
      

#########################
#      Main Driver      #
#########################
    def run(self):
        self.readOpt()
        self.sanityCheck()
        self.readConformers()
        self.pruneByEnergy()
        self.align()
        self.plot()
         
if __name__ == '__main__':
    p = plotConformers()
    p.run()
