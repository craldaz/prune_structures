from scipy.spatial import distance
from scipy.optimize import minimize
import numpy as np
import getopt
import pybel
import sys
import subprocess
import math
from collections import namedtuple 

Torsion=namedtuple('Torsion',['IDs','angles'])

class genConformer:
    def __init__(self, MolName="", TorsionFileName="", OutputDir="",createBond=[],Smi=""):
        self.molName_ = MolName 
        self.torsionFileName_ = TorsionFileName
        self.conformer_ =0
        self.maxConformer_ = 1
        self.currentTorsion_ = []
        self.currentTargetTorsion_ = []
        self.torsionList_ = []
        self.indexed_ = 1
        self.myOBMol_ = None
        self.refSmi_ =Smi
        #IO names
        self.outDir_ = OutputDir
        self.baseNames_ = []
        self.xyzNames_ = []
        #Things related to formatting
        self.spacePerAtom_ = 0
        self.spacePerTorsion_ = 0
        self.spacePerConformer_ = 0
        # Parameters
        self.FailedTorsionThre_ = 1.0
        # pre-rotation action
        self.createBond_=createBond

    def usage(self):
        print "Usage: python rotate.py -i [xyz_file] -t [torsion_angle_specification_file]"
        print ''' 
                  -i   xyz file (Required)
                  -t   torsion angle specification (Required)
                       Format: 
                       atom1 atom2 atom3 atom4  angle1 angle2 angle3 angle4 ....
                       '''

    def readOpt(self):
        try:
          optlist,args = getopt.getopt(sys.argv[1:],'i:t:')
        except getopt.GetoptError as err:
          print str(err)
          self.usage()
          sys.exit(2)   
        for o,a in optlist:
             if o == "-i":
               self.molName_ = a
             elif o == "-t":
               self.torsionFileName_ = a
             else:
               assert False, "unhandled option"
    def readTorsionFile(self):
        data=open(self.torsionFileName_).readlines()
        try:
           self.indexed_ = int(data[0].strip('\n'))
        except:
            print 'Wrong format on 1st line of Torsion file. Must be 0 or 1'
            sys.exit(2)
        for line in data[1:]:
            items = line.strip('\n').split()
            try:
                IDs=[int(items[i]) for i in range(0,4)]
            except IndexError as err:
                print str(err)
                print 'Line ended unexpectedly when reading torsion IDs.'
                sys.exit(2)
            except ValueError:
                print 'First 4 columns of torsion file must be integers.'
                sys.exit(2)
             
            try:
                angles=[ float(item) for item in items[4:] ]
            except IndexError as err:
                print str(err)
                print 'Line ended unexpectedly when reading torsion angle values.'
                sys.exit(2)
            except ValueError:
                print 'Columns 4 to last of torsion file must be floats.'
                sys.exit(2)
            self.torsionList_.append(Torsion(IDs,angles))

        print 'Original atom IDs are '+str(self.indexed_) +'-indexed'
            

    def findTorsions(self):
        self.myOBMol_.FindTorsions()
        rotorList=[]
        for torsion in pybel.ob.OBMolTorsionIter(self.myOBMol_): # Notice that the torsion indices are 0-indexed
            if not (self.myOBMol_.GetAtom(torsion[1]+1).IsInRing() and self.myOBMol_.GetAtom(torsion[2]+1).IsInRing()):
                 duplicated=False
                 for rotor in rotorList:
                     if rotor[0][1] == torsion[1] and rotor[0][2] == torsion[2]:
                         duplicated=True
                         print 'rotor duplicated',rotor
                         break
                 if duplicated == False:
                     mybond=self.myOBMol_.GetBond(torsion[1]+1,torsion[2]+1)
                     mybondOrder= mybond.GetBondOrder()
                     rotorList.append([torsion,mybondOrder])
        for rotor in rotorList:
            [torsion,bondOrder]=rotor
            IDs=[ torsion[i]+1 for i in range(0,4) ]
            if bondOrder == 1:
               angles=[-120, 0, 120]
            elif bondOrder == 2:
               angles=[0, 180]
            self.torsionList_.append(Torsion(IDs,angles))
        print self.torsionList_
               
        

    def sanityCheck(self, verbose=False):
        if self.molName_ =="":
              assert False, "xyz file must be provided!"
        else:
              if verbose == True:
                print "xyz file:",self.molName_
              mymol=pybel.readfile('xyz',self.molName_).next()
              self.myOBMol_ = mymol.OBMol

        if self.torsionFileName_ =="":
              #assert False, "torsion file must be provided!"
              self.findTorsions()
        else:
              if verbose == True: 
                print "torsion file:",self.torsionFileName_
              self.readTorsionFile()

        #Check whether the torsions specified are legal based on the xyz structure
        for torsion in self.torsionList_:
            for ID in torsion.IDs:
                if ID < 0 or ID > self.myOBMol_.NumAtoms():
                    assert False, 'The torsion specified is not valid because atom index out of list'+str(IDs)
                    break

         #Initialize class member values
        self.spacePerAtom_ = int(np.ceil(np.log10(self.myOBMol_.NumAtoms())))
        self.spacePerTorsion_ = 4*(self.spacePerAtom_+1)-1
        self.currentTorsion_ = [0]*len(self.torsionList_)
        self.currentTargetTorsion_ = [0]*len(self.torsionList_)
        for torsion in self.torsionList_:
            self.maxConformer_ *= len(torsion.angles)
        self.spacePerConformer_ = int(np.ceil(np.log10(self.maxConformer_)))
        self.separator_ = ' '.rjust(self.spacePerAtom_)
        print 'Maximum number of conformers going to be generated: ', self.maxConformer_

    def preRotateAction(self):
        if len(self.createBond_)>0:
           for bond in self.createBond_:
               [beginIdx,endIdx,order]=bond
               if self.indexed_ == 0:
                 status=self.myOBMol_.AddBond(beginIdx+1,endIdx+1,order)
                 print 'Added bond (0-indexed) ',beginIdx,'=>',endIdx,'order',order,'success',status
               else:
                 status=self.myOBMol_.AddBond(beginIdx,endIdx,order)
                 print 'Added bond (1-indexed) ',beginIdx,'=>',endIdx,'order',order,'success',status

    #Return the difference considering the periodicity of degrees
    def getAngleDifference(self, angle1,angle2):
        difference=[math.fabs(angle1-angle2), math.fabs(angle1+360-angle2), math.fabs(angle1-360-angle2)]
        return min(difference)
                
    def checkGeometry(self,mymol):
        difference=[ self.getAngleDifference(self.currentTorsion_[i] , self.currentTargetTorsion_[i]) for i in range(0,len(self.torsionList_))]
        # 1. Check whether the torsion angle cannot be reached
        if max(difference) > self.FailedTorsionThre_:
             torsionLevel= difference.index(max(difference))
             torsion = self.torsionList_[torsionLevel]
             angle = self.currentTargetTorsion_[torsionLevel]
             print 'Current torsion invalid because the target angle ',angle, 'of the torsion ',torsion.IDs,' cannot be reached'
             print 'Skip this conformer'
             return False
        # 2. Check whether the system is collisions between atoms
        if self.refSmi_ !="":
           refMol=pybel.readstring("smi",self.refSmi_)
           refFp =refMol.calcfp()
           outname='temp.xyz'
           mymol.write('xyz',outname,overwrite=True)
           mol = pybel.readfile("xyz",'temp.xyz').next()
           subprocess.call('rm temp.xyz',shell=True)
           fp = mol.calcfp()
           TanimotoCoeff = refFp|fp
           if  TanimotoCoeff < 1:
             print 'Current torsion invalid because current geometry has different chemical structure than the original one (caused by some collisions between atoms) '
             print "TanimotoCoeff:",TanimotoCoeff
             print 'Skip this conformer'
             return False
        return True
           

    def printCurrentTorsion(self):
        info='conformer '+str(self.conformer_).rjust(self.spacePerConformer_)
        for angle in self.currentTorsion_:
           info+="{0:.1f}".format(angle).rjust(self.spacePerTorsion_)+self.separator_
        info+='\n'
        sys.stdout.write(info)
        self.myOBMol_.SetTitle('0')

    def singleTorsion(self,torsionLevel):
        if torsionLevel == len(self.torsionList_):
              self.printCurrentTorsion()
              mymol=pybel.Molecule(self.myOBMol_)
              if self.checkGeometry(mymol):
                 #basename='c'+str(self.conformer_).zfill(self.spacePerConformer_)
                 basename='conf'+str(self.conformer_)
                 fullbasename=self.outDir_+basename
                 outname=fullbasename+'.xyz'
                 self.baseNames_.append(basename)
                 self.xyzNames_.append(outname)
                 mymol.write('xyz',outname,overwrite=True)
                 self.conformer_+=1
              return
        else:
              torsion=self.torsionList_[torsionLevel]
              IDs=torsion.IDs
              angles=torsion.angles
              atoms =[ self.myOBMol_.GetAtom(IDs[i]) for i in range(0,4)]
              for angle in angles:
                 targetAngle=math.radians(angle)
                 self.myOBMol_.SetTorsion(atoms[0],atoms[1],atoms[2],atoms[3],targetAngle)
                 realAngle= self.myOBMol_.GetTorsion(IDs[0],IDs[1],IDs[2],IDs[3])
                 self.currentTorsion_[torsionLevel] = realAngle
                 self.currentTargetTorsion_[torsionLevel] = angle
                 self.singleTorsion(torsionLevel+1)

    def printTitle(self):
        print 'List of torsion angles for each conformer:'
        if(self.indexed_ == 1):
          sys.stdout.write('(1-indexed)'.rjust(10+ self.spacePerConformer_))
        else:
          sys.stdout.write('(0-indexed)'.rjust(10+ self.spacePerConformer_))
        for torsion in self.torsionList_:
            IDs=torsion.IDs
            for i in IDs[:-1]:
               sys.stdout.write(str(i).rjust(self.spacePerAtom_)+'-')
            sys.stdout.write(str(IDs[-1]).rjust(self.spacePerAtom_)+self.separator_)
        sys.stdout.write('\n')

    def adjustAtomIDs(self):
        if(self.indexed_ == 0):
           print 'Adjust atom IDs to 1-indexed to be further processed'
           for torsion in self.torsionList_:
               for i in range(0,4):
                  torsion.IDs[i] += 1
           self.indexed_ = 1
           self.printTitle()

    def summary(self):
        print 'Finished processing all torsion angles'
        print 'Ideal number of conformers to be generated: ', self.maxConformer_
        print 'Number of conformers actually generated: ', self.conformer_
      

    def rotate(self):
        self.printTitle()
        self.adjustAtomIDs()
        self.singleTorsion( 0)

    def run(self):
        self.readOpt()
        self.sanityCheck(verbose=True)
        self.preRotateAction()
        self.rotate()
        self.summary()

    def runAsModule(self):
        self.sanityCheck(verbose=True)
        self.preRotateAction()
        self.rotate()
        self.summary()
   
    def getNumAtoms(self):
        return self.myOBMol_.NumAtoms()
    
    def getXyzNames(self):
        return self.xyzNames_
   
    def getBaseNames(self):
        return self.baseNames_
if __name__ == '__main__':
    me=genConformer()
    me.run()
       
