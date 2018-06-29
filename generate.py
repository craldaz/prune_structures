import subprocess
import sys
import os
from time import sleep
import getopt
import pybel
import numpy as np
import math
import pickle
import rotate
import itertools

# Function: printLine
def printLine():
    print "------------------------------------------------------"
    print '\n'

# Function: com2xyz
"""------------------
 Transform the Gaussian com input file 
 intput xyz file that can be read by TeraChem
"""
def com2xyz(infile):
    datain=open(infile).readlines()
    outfilename=str(infile).replace('com','xyz')
    #First determine the number of atoms in the mol and where the coords start in inputdata
    startline=-1
    endline=-1
    for i in range(0,len(datain)):
      line=datain[i].strip(" ").split()
      if len(line)>0 and line[0]=="conformer":
          startline=i
          break
    if startline<0:
      print "Wrong input file format! lacking start line"
      exit
    for i in range(startline,len(datain)):
      line=datain[i].strip(" ").split()
      if len(line)>0 and line[0]=="Extbasis":
          endline=i
          break
    if endline<0:
      print "Wrong input file format! lacking endline"
      exit
    #startline is the line# of the 1st coordinate line. Endline is the line# of the last coord line
    startline+=3
    endline-=2
    natom = endline+1-startline
    outfile=open(outfilename,"w")
    outfile.write(str(natom)+"\n\n")
    for i in range(startline,endline+1):
       outfile.write(datain[i])
    outfile.close()

# Function: jobFinished
''' ------------------------------
 Check whether all the minimization jobs are finished
 This can be run in several mode
 mode=1: check whether the minimization geometry is generated
 mode=2: check whether the job number no longer exist in the job list
 Mode 2 should be safer, as sometimes one job fails to generate the
 wanted file and this function in Mode 1 will return false forever
'''
def jobFinished(jobIDs,username):
    status=True
    proc=subprocess.Popen(['squeue','-u',username],
                           stdout=subprocess.PIPE)
    runningJobs=proc.communicate()[0]
    for myJobID in jobIDs:
        if runningJobs.find(myJobID)>=0:
           status=False
           break
        
    return status

#Function readXyzEnergy(xyzfile)
def readXyzEnergy(xyzfile):
    data=open(xyzfile).readlines()
    line=data[1].split()
    if len(line) < 1:
       print 'Wrong XYZ file format!'
       sys.exit(1)
    return float(line[0])

class confgen:
    def __init__(self):
        self.molName_ = ""
        self.torsionFileName_ = ""
        self.generator_ = None
        self.workDir_ = ""
        self.userName_ = ""
        self.checkStatusTime_ = 5
        self.restart_ = False
        self.verbose_ = False
        self.step_ = 0
        self.xyzNames_ = []
        self.xyzDir_ = ""
        self.baseNames_ = []
        self.jobIDs_ = []
        self.notConvergedNames_ = []
        self.createBond_ = []
        self.indexed_ = 0
        self.numAtoms_ = 0
        self.modFile_ = ""
        self.no_sorting_ = False
        self.no_opt_ = False
        self.refSmi_ =""
        self.localopt_ = False
    
    def printHeader(self):
        print """ 
        ##########################
        # Starting AUTO NEB run  #
        ########################## """       
    
    def usage(self):
        print "Usage: python generate.py -i [xyz_file] -t [torsion_angle_specification_file] -m [xyz modification file]  [-r] [-v] [-T [check_time]] -d [unoptimized_xyz_file_path] --no_sorting --no_opt --smi [reference smi]"
        print ''' 
                  -i   xyz file (Required)
                  -t   torsion angle specification (Required)
                       Format: 
                       atom1 atom2 atom3 atom4  angle1 angle2 angle3 angle4 ....
                  -r   Restart mode
                  -v   Verbose output
                  -m   xyz modification file
                  -T   Time interval to check whether the minimization jobs 
                       are finished
                  -b   Bonds to be created
                  --smi reference SMILE string to prune wrong molecules
                  --no_sorting: do no sort the minimized structure
                  --no_opt: do not optimize the structures (usually used when we just want to create some new structures)
                  --localopt: do local optimization with pybel's internal force field after doing modification to the structures (currently only works with mode=3) '''

    def readOpt(self):
        print 'Argument List:', str(sys.argv)
        try:
          optlist,args = getopt.getopt(sys.argv[1:],'vri:t:T:d:b:m:',['no_sorting','no_opt','smi=','localopt'])
        except getopt.GetoptError as err:
          print str(err)
          self.usage()
          sys.exit(2)   
        for o,a in optlist:
             if o == "-i":
               self.molName_ = a
             elif o == "-t":
               self.torsionFileName_ = a
             elif o == "-r":
               self.restart_ = True
             elif o == "-v":
               self.verbose_ = True
             elif o == "-T":
               self.checkStatusTime_ = float(a)
             elif o == "-d":
               self.xyzDir_ = a
             elif o == "-m":
               self.modFile_ = a
             elif o == "-b":
               bonds= a.split(',')
               print len(bonds)
               for i in range(0,len(bonds)/3):
                  self.createBond_.append([int(bonds[3*i]),int(bonds[3*i+1]),int(bonds[3*i+2])])
               print self.createBond_
             elif o == "--no_sorting":
                 self.no_sorting_ =True
             elif o == "--no_opt":
                 self.no_opt_ =True
             elif o == "--localopt":
                 self.localopt_ =True
             elif o == "--smi":
               self.refSmi_ = a
             else:
               assert False, "unhandled option"

    def printSanityCheck(self,term,val,append=""):
        print term.ljust(40),val,append
     
    def sanityCheck(self):
        self.printHeader()
        self.workDir_ = subprocess.Popen(['pwd'],
                         stdout=subprocess.PIPE).communicate()[0].strip('\n')
        self.userName_ = subprocess.Popen(['whoami'],
                         stdout=subprocess.PIPE).communicate()[0].strip('\n')
        if self.restart_ == False and self.modFile_=="":
           print "Running AUTO NEB from scratch."
           if self.molName_ =="":
                 assert False, "xyz file must be provided!"
           else:
                 self.printSanityCheck("xyz file: ",self.molName_)
           
           if self.torsionFileName_ =="":
                 #assert False, "torsion file must be provided!"
                 print "No torsion file provided. Torsions will be generated automatically."
           else:
                 self.printSanityCheck("torsion file:",self.torsionFileName_)
           if self.xyzDir_ == "":
                 self.xyzDir_ = "./xyz/"
        elif self.modFile_!="": 
           print "Running AUTO NEB by modifying existing xyz." 
           if self.xyzDir_ == "":
                 assert False, "xyz Directory must be provided!"
           else:
                 self.printSanityCheck("xyz directory: ",self.xyzDir_)
        else:
           print "Running AUTO NEB in RESTART MODE." 
           try:
              input=open('autoneb.chk','r')
           except IOError as e:
              print 'Restart checkpoint file autoneb.chk does not exist'
              print "I/O error({0}): {1}".format(e.errno, e.strerror)
           self.notConvergedNames_=pickle.load(input)
           print "Successfully loaded the checkpoint file for restarting minimization"

       
        self.printSanityCheck('My work directory is ',self.workDir_)
        self.printSanityCheck('My username is ',self.userName_)
        self.printSanityCheck('Generated conformers (before optimization) will be written to ', self.xyzDir_)
        self.printSanityCheck('Job status will be checked every',self.checkStatusTime_,'seconds')
        if self.no_sorting_ == True:
           self.printSanityCheck('Minimized structures will not be sorted','')
        if self.no_opt_ == True:
           self.printSanityCheck('Generated conformers will not be optimized','')
        if self.refSmi_ !="":
           self.printSanityCheck('Generated conformers will be checked against SMILE string:',self.refSmi_)
        printLine()
        self.step_+=1

#######################################
# 1) Generate conformers from scratch #
#######################################
    def generateXyz(self):
        print self.step_,'. Rotating torsion angles to generate conformers...'
        cmd = 'mkdir -p '+ self.xyzDir_
        subprocess.call(cmd,  shell=True)
        self.generator_ = rotate.genConformer(self.molName_, self.torsionFileName_, self.xyzDir_, self.createBond_,self.refSmi_)
        self.generator_.runAsModule()
        self.baseNames_ = self.generator_.getBaseNames()
        self.xyzNames_ = self.generator_.getXyzNames()
        print 'Finished generating conformers!'
        printLine()
        self.step_ +=1

#############################################
# 2) Read conformers and modify structures  # 
#    (for modification mode)                #
#############################################
    def getFileNames(self,fileDirectory):
        cmd=['ls',fileDirectory]
        proc=subprocess.Popen(cmd,  
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE
                              )
        stdout_value = proc.communicate()[0] #Get the stdout output
        names = str(stdout_value).split()
        names = [ name.strip('.xyz') for name in names ]
        return names

    def readConformers(self):
        self.baseNames_ =self.getFileNames(self.xyzDir_)
        numXyz=len(self.baseNames_)
        if numXyz < 1:
           print 'No xyz file found in current directory. Nothing to do.'
           return
        
        print 'Number of conformers read in:',numXyz
        for name in self.baseNames_:
            xyzName = self.xyzDir_ +'/'+name+'.xyz'
            self.xyzNames_.append(xyzName)

    def calcCoordsSp2(self, myMol, i0, i1, i2, bondlength):
        hcoordsV3i = pybel.ob.vector3(0,0,0)
        success = False
        coords = [ iatom.coords for iatom in myMol.atoms ]
        b1=np.array(coords[i0])-np.array(coords[i1])
        b2=np.array(coords[i0])-np.array(coords[i2])
        v= b1 +b2
        vnorm = np.linalg.norm(v)
        if vnorm == 0:
            print "Error in calculating position of Hydrogen"
        else:
            v /= vnorm
            hcoords = np.array(coords[i0]) + bondlength*v
            hcoordsV3 = pybel.ob.vector3(hcoords[0],hcoords[1],hcoords[2])
            success = True
        return [hcoordsV3,success]

    def calcCoordsOH(self, myMol, i0, i1, i2, bondlength):
        hcoordsV3i = pybel.ob.vector3(0,0,0)
        success = False
        coords = [ iatom.coords for iatom in myMol.atoms ]
        v=np.array(coords[i1])-np.array(coords[i2])
        vnorm = np.linalg.norm(v)
        if vnorm == 0:
            print "Error in calculating position of Hydrogen"
        else:
            v /= vnorm
            hcoords = np.array(coords[i0]) + bondlength*v
            hcoordsV3 = pybel.ob.vector3(hcoords[0],hcoords[1],hcoords[2])
            success = True
        return [hcoordsV3,success]

    def addHydrogenSp2(self, myMol, HPosition):
        status=False
        [Htype,i0,i1,i2,i3,bondlength]=HPosition
        [hcoordsV3,success] = self.calcCoordsSp2(myMol, i0, i1, i2, bondlength)
        if success == False:
           status=False
        else:
            newOBAtom = pybel.ob.OBAtom()
            newOBAtom.SetVector(hcoordsV3)
            newOBAtom.SetAtomicNum(1)
            status = myMol.OBMol.AddAtom(newOBAtom)
        return status

    def moveHydrogenSp2(self, myMol, HPosition):
        [Htype,iRemove,i0,i1,i2,i3,bondlength]=HPosition
        #Get the hydrogen at iRemove
        atomToMove = myMol.OBMol.GetAtom(iRemove+1)
        #calculate the new position of the hydrogen (attach to i0)
        [hcoordsV3,success] = self.calcCoordsSp2(myMol, i0, i1, i2, bondlength)
        if success == False:
           status=False
        else:
           status=atomToMove.SetVector(hcoordsV3)
        return status

    def moveHydrogenOH(self, myMol, HPosition):
        [Htype,iRemove,i0,i1,i2,i3,bondlength]=HPosition
        #Get the hydrogen at iRemove
        atomToMove = myMol.OBMol.GetAtom(iRemove+1)
        #calculate the new position of the hydrogen (attach to i0)
        [hcoordsV3,success] = self.calcCoordsOH(myMol, i0, i1, i2, bondlength)
        if success == False:
           status=False
        else:
           status=atomToMove.SetVector(hcoordsV3)
        return status
    #Add H atom to specified positions to a single xyz file
    # Atoms are 0-indexed
    def addHydrogenSingleMol(self, xyzName, HPositionList):
        myMol = pybel.readfile('xyz', xyzName).next()
        statusTotal = True
        for HPosition in HPositionList:
           #Calculate new hydrogen position and add hydrogen
           Htype=HPosition[0]
           if Htype=="sp2":
              status = self.addHydrogenSp2(myMol, HPosition)
              #Print out Hydrogen addition result
              if status == False:
                 if self.indexed_ == 0:
                   print 'Failed to add sp2 H atom to molecule: ',xyzName ,', at atom ',i0,'(',i1,i2,'), bondlength:',round(bondlength,2)
                 else:
                   print 'Failed to add sp2 H atom to molecule: ',xyzName ,', at atom ',i0+1,'(',i1+1,i2+1,'), bondlength:',round(bondlength,2)
                   statusTotal = False
                   break
           else:
              #TODO: implement this part
              print "sp3 H adding not implemented yet"
              statusTotal = False
              break
        # delete the original xyz file
        subprocess.call(['rm','-r',xyzName])
        # if succeeded, overwrite and get new xyz
        if statusTotal == True:
           output=pybel.Outputfile("xyz",xyzName,overwrite=True)
           output.write(myMol)
           output.close()
        return statusTotal

    #Move H atom from a position to specified position to a single xyz file
    # Atoms are 0-indexed
    def moveHydrogenSingleMol(self, xyzName, HPositionList):
        myMol = pybel.readfile('xyz', xyzName).next()
        statusTotal = True
        for HPosition in HPositionList:
           #Calculate new hydrogen position and add hydrogen
           Htype=HPosition[0]
           if Htype=="sp2":
             status = self.moveHydrogenSp2(myMol, HPosition)
             #Print out Hydrogen move result
             if status == False:
                if self.indexed_ == 0:
                  print 'Failed to move H atom from',iRemove,' to',i0,' in molecule: ',xyzName, ' bondlength:',round(bondlength,2)
                else:
                  print 'Failed to move H atom from',iRemove+1,' to',i0+1,' in molecule: ',xyzName ,' bondlength:',round(bondlength,2)
                statusTotal =False
                break
           elif Htype=="oh":
             status = self.moveHydrogenOH(myMol, HPosition)
             if status == False:
                if self.indexed_ == 0:
                  print 'Failed to move H atom from',iRemove,' to',i0,' in molecule: ',xyzName, ' bondlength:',round(bondlength,2)
                else:
                  print 'Failed to move H atom from',iRemove+1,' to',i0+1,' in molecule: ',xyzName ,' bondlength:',round(bondlength,2)
                statusTotal =False
                break
           else:
             #TODO: implement this part
             print "sp3 H adding not implemented yet"
             statusTotal =False
             break
        # delete the original xyz file
        subprocess.call(['rm','-r',xyzName])
        # if succeeded, overwrite and get new xyz
        if statusTotal == True:
           output=pybel.Outputfile("xyz",xyzName)
           output.write(myMol)
           output.close()
        return statusTotal
      
    def addHydrogens(self):
        inp = open(self.modFile_).readlines()
        HPositionList=[]
        for i in range(1,len(inp)):
            line=inp[i]
            e=line.split()
            [Htype,i0,i1,i2,i3,bondlength]=["",-1,-1,-1,-1,-1]
            if e[0] == "sp2":
                [Htype,i0,i1,i2,bondlength]=["sp2",int(e[1]),int(e[2]),int(e[3]),float(e[4])]
            elif e[0]=="sp3":
                [Htype,i0,i1,i2,i3,bondlength]=["sp3",int(e[1]),int(e[2]),int(e[3]),int(e[4]),float(e[5])]
            else:
                print "Unexpected format in xyz modification file: ",self.modFile_,", line: ",i+1
                sys.exit(2)
            # adjust indices if 1-indexed
            if self.indexed_ == 1:
                i0-=1
                i1-=1
                i2-=1
                if i3 > -1:
                     i3-=1
            HPositionList.append([Htype,i0,i1,i2,i3,bondlength])
        statusTotal=True
        for xyzName in self.xyzNames_:
            status=self.addHydrogenSingleMol(xyzName,HPositionList)
            if status == True:
               print "Successfully added hydrogens to ",xyzName
            else:
               print "Failed to add hydrogens to ",xyzName
            statusTotal = statusTotal and status
        return statusTotal

    def moveHydrogens(self):
        inp = open(self.modFile_).readlines()
        HPositionList=[]
        for i in range(1,len(inp)):
            line=inp[i]
            e=line.split()
            [Htype,iRemove,i0,i1,i2,i3,bondlength]=["",-1,-1,-1,-1,-1,-1]
            if e[0] == "sp2":
                [Htype,iRemove,i0,i1,i2,bondlength]=["sp2",int(e[1]),int(e[2]),int(e[3]), int(e[4]),float(e[5])]
            elif e[0]=="sp3":
                [Htype,iRemove,i0,i1,i2,i3,bondlength]=["sp3",int(e[1]),int(e[2]),int(e[3]),int(e[4]),int(e[5]),float(e[6])]
            elif e[0]=="oh":
                [Htype,iRemove,i0,i1,i2,bondlength]=["oh",int(e[1]),int(e[2]),int(e[3]), int(e[4]),float(e[5])]
            else:
                print "Unexpected format in xyz modification file: ",self.modFile_,", line: ",i+1
                sys.exit(2)
            # adjust indices if 1-indexed
            if self.indexed_ == 1:
                iRemove-=1
                i0-=1
                i1-=1
                i2-=1
                if i3 > -1:
                     i3-=1
            HPositionList.append([Htype,iRemove,i0,i1,i2,i3,bondlength])
        statusTotal=True
        for xyzName in self.xyzNames_:
            status = self.moveHydrogenSingleMol(xyzName,HPositionList)
            if status == True:
               print "Successfully moved hydrogens in ",xyzName
            else:
               print "Failed to move hydrogens in ",xyzName
            statusTotal = statusTotal and status
        return statusTotal

    # i0, i1 are 0-indexed atom indices
    def modifyBondsSingleMol(self, xyzName,bondList):
        myMol = pybel.readfile('xyz', xyzName).next()
        statusTotal = True
        for bond in bondList:
            [bondLength,i0,i1]=bond
            bondToChange = myMol.OBMol.GetBond(i0+1,i1+1)
            if bondToChange == None:
               print "WARNING! The specified bond ",i0,"-",i1," does not exist"
               print "Creating bond ",i0,"-",i1," and then process the bondlength modification request"
               bondOrder =1
               addBondStatus=myMol.OBMol.AddBond(i0+1,i1+1,bondOrder)
               if addBondStatus == False:
                  print "Failed to add specified bond ",i0,"-",i1
                  statusTotal = False
                  break
               else:
                  bondToChange = myMol.OBMol.GetBond(i0+1,i1+1)
            atomToFix= myMol.OBMol.GetAtom(i0+1)
            bondToChange.SetLength(atomToFix,bondLength)
        if self.localopt_ == True:
           myMol.localopt(steps=500) 
         
        # delete the original xyz file
        subprocess.call(['rm','-r',xyzName])
        # if succeeded, overwrite and get new xyz
        if statusTotal == True:
           output=pybel.Outputfile("xyz",xyzName)
           output.write(myMol)
           output.close()
        return statusTotal

#d1,d2 are 0-indexed atom numbers
    def moveGroupSingleMol(self, xyzName,modList):
        myMol = pybel.readfile('xyz', xyzName).next()
        statusTotal = True
        coords = [ iatom.coords for iatom in myMol.atoms ]
        for mod in modList:
            [dist,d1,d2,groupAtoms]=mod
            if len(groupAtoms) <1:
                 print "Nothing to move along direction ",d1,"-->",d2
            [dist,d1,d2,groupAtoms]=mod
            v1=np.array(coords[d1])
            v2=np.array(coords[d2])
            v= v2-v1
            vnorm=np.linalg.norm(v)
            v/=vnorm
            for i in groupAtoms:
                newcoord = np.array(coords[i]) + dist*v
                newcoordV3 = pybel.ob.vector3(newcoord[0],newcoord[1],newcoord[2])
                atomToMove = myMol.OBMol.GetAtom(i+1) 
                atomToMove.SetVector(newcoordV3)
        # delete the original xyz file
        subprocess.call(['rm','-r',xyzName])
        # if succeeded, overwrite and get new xyz
        if statusTotal == True:
           output=pybel.Outputfile("xyz",xyzName)
           output.write(myMol)
           output.close()
        return statusTotal

    def moveGroup(self):
        inp = open(self.modFile_).readlines()
        modList=[]
        for i in range(1,len(inp)):
            line=inp[i]
            e=line.split()
            [dist,d1,d2]=[0.0,-1,-1]
            groupAtoms=[]
            try:
                [dist,d1,d2]=[float(e[0]),int(e[1]),int(e[2])]
                for item in e[3:]:
                    groupAtoms.append(int(item))
            except ValueError:
                print "Unexpected format in xyz modification file: ",self.modFile_,", line: ",i+1
                sys.exit(2)
            # adjust indices if 1-indexed
            if self.indexed_ == 1:
                d1-=1
                d2-=1
                for i in range(0,len(groupAtoms)):
                    groupAtoms[i]-=1
            modList.append([dist,d1,d2,groupAtoms])
        statusTotal=True
        for xyzName in self.xyzNames_:
            status = self.moveGroupSingleMol(xyzName,modList)
            if status == True:
               print "Successfully moved atoms in ",xyzName
            else:
               print "Failed to move atoms in ",xyzName
            statusTotal = statusTotal and status
        return statusTotal
   
    def rotateAroundAxis(self, oldcoord, vcenter, axis, theta):
        #calculate the rotation matrix
        alpha = np.deg2rad(theta)
        cosA = np.cos(alpha)
        sinA = np.sin(alpha)
        [x,y,z]=axis
        M=np.array([[   cosA      + (1-cosA)*x*x, (1-cosA)*x*y-   sinA   *z, (1-cosA)*x*z+      sinA*y],
                    [(1-cosA)*y*x +    sinA   *z,    cosA     +(1-cosA)*y*y, (1-cosA)*y*z-      sinA*x],
                    [(1-cosA)*z*x -    sinA   *y, (1-cosA)*z*y+    sinA  *x,    cosA     +(1-cosA)*z*z]])
        ca = oldcoord-vcenter
        da = ca- np.dot(ca,axis)
        d =  oldcoord - da
        db = np.dot(M,da)
        b = d + db 
        print "Rotation around axis", axis, " centered at ", vcenter,": ", oldcoord, "-->", b
        return b


    def rotateGroupSingleMol(self, xyzName,modList):
        myMol = pybel.readfile('xyz', xyzName).next()
        statusTotal = True
        coords = [ iatom.coords for iatom in myMol.atoms ]
        [angle, groupA, groupB] = modList
        if len(groupA) <1:
             print "Nothing to rotate along direction "
             return statusTotal
        #Get the center of set B
        vcenter = np.array([0.0,0.0,0.0])
        for i in groupB:
           vcenter+=np.array(coords[i])
        vcenter/=len(groupB)
        #Get the axis of the plan defined by set B
        v1 = np.array(coords[groupB[0]]) - np.array(coords[groupB[1]])
        v2 = np.array(coords[groupB[2]]) - np.array(coords[groupB[1]])
        axis = np.cross(v1,v2)
        vnorm=np.linalg.norm(axis)
        axis/=vnorm
        #Now do the rotation
        for i in groupA:
            oldcoord = np.array(coords[i])
            newcoord = self.rotateAroundAxis(oldcoord, vcenter, axis, angle)
            newcoordV3 = pybel.ob.vector3(newcoord[0],newcoord[1],newcoord[2])
            atomToMove = myMol.OBMol.GetAtom(i+1) 
            atomToMove.SetVector(newcoordV3)
        # delete the original xyz file
        subprocess.call(['rm','-r',xyzName])
        # if succeeded, overwrite and get new xyz
        if statusTotal == True:
           output=pybel.Outputfile("xyz",xyzName)
           output.write(myMol)
           output.close()
        return statusTotal

    def rotateGroup(self):
        inp = open(self.modFile_).readlines()
        if(len(inp))<4:
           print "Wrong format of xyz modification file for mode 8: rotate a group. The file should have 4 lines"
           sys.exit(2)
        elif(len(inp))>4:
           print "Contents after the 4th line of the xyz modification file for mode 8 will be ignored"

        angle=float(inp[1].strip('\n').split()[0])
        line = inp[2].strip('\n').split()
        groupA = []
        groupB = []
        try:
            groupA = [ int(x) for x in line ]
        except ValueError:
            print "Unexpected format in xyz modification file: ",self.modFile_,", line: 3"
        line = inp[3].strip('\n').split()
        try:
            groupB = [ int(x) for x in line ]
        except ValueError:
            print "Unexpected format in xyz modification file: ",self.modFile_,", line: 4"
            sys.exit(2)
        if len(groupB) <3:
           print "For modification mode 8, at least 3 atoms are needed for set B"
        # adjust indices if 1-indexed
        if self.indexed_ == 1:
            for i in range(0,len(groupA)):
                groupA[i]-=1
            for i in range(0,len(groupB)):
                groupB[i]-=1
        modList=[ angle,groupA,groupB]
        print "Modification mode 8: rotate atoms around a axis"
        print "Rotation degree: ", angle
        print "Rotation atoms: " , groupA
        print "Rotation axis is the plane axis of atoms: ", groupB
        statusTotal=True
        for xyzName in self.xyzNames_:
            status = self.rotateGroupSingleMol(xyzName,modList)
            if status == True:
               print "Successfully roated atoms in ",xyzName
            else:
               print "Failed to rotated atoms in ",xyzName
            statusTotal = statusTotal and status
        return statusTotal

    def modifyBonds(self):
        inp = open(self.modFile_).readlines()
        bondList=[]
        for i in range(1,len(inp)):
            line=inp[i]
            e=line.split()
            [bondLength,i0,i1]=[-1,-1,-1]
            try:
                [bondLength,i0,i1] = [float(e[0]),int(e[1]),int(e[2])]
            except ValueError:
                print "Unexpected format in xyz modification file:",self.modFile_,",line:",i+1
                sys.exit(2)
            # adjust indices if 1-indexed
            if self.indexed_ == 1:
                i0-=1
                i1-=1
            bondList.append([bondLength, i0, i1])
        statusTotal=True
        for xyzName in self.xyzNames_:
            status = self.modifyBondsSingleMol(xyzName,bondList)
            if status == True:
               print "Successfully modified bonds in ",xyzName
            else:
               print "Failed to modify bonds in ",xyzName
            statusTotal = statusTotal and status
        return statusTotal

    def removeAtoms(self):
        inp = open(self.modFile_).readlines()
        atomIDList=[]
        line=inp[1]
        e=line.strip('\n').split()
        for entry in e:
          try:
              atomIDList.append(int(entry))
          except ValueError:
              print "Unexpected format in xyz modification file:",self.modFile_,",line:2"
              sys.exit(2)
        # adjust indices if 0-indexed
        if self.indexed_ == 0:
           for i in range(len(atomIDList)):
                atomIDList[i]+=1 
        statusTotal=True
        for xyzName in self.xyzNames_:
            myMol = pybel.readfile('xyz', xyzName).next()
            atomList = []
            myMol.OBMol.BeginModify()
            for atom in pybel.ob.OBMolAtomIter(myMol.OBMol):
                if atom.GetIdx() in atomIDList:
                   atomList.append(atom)
            for atom in atomList:
                   status = myMol.OBMol.DeleteAtom(atom)
                   if status ==False:
                      print "Failed to delete atom in ",xyzName
                      statusTotal = False
                      break
            myMol.OBMol.EndModify()
            if statusTotal == False:
               break
            else:
               print "successfully deleted atoms in ",xyzName
               output=pybel.Outputfile("xyz",xyzName,overwrite=True)
               output.write(myMol)
               output.close()
               
        return statusTotal

    def dissociateC2H4(self, moveAway):
        status=True
        inp = open(self.modFile_).readlines()
        line=inp[1]
        e=line.split()
        [bondLength,iO,iH,iC1,iH11,iH12,iH13,iC2,iH21,iH22]=[-1.0,-1,-1,-1,-1,-1,-1,-1,-1,-1]
        try:
            [bondLength,iO,iH,iC1,iH11,iH12,iH13,iC2,iH21,iH22] = [float(e[0]),int(e[1]),int(e[2]), int(e[3]),int(e[4]),int(e[5]),int(e[6]),int(e[7]), int(e[8]), int(e[9])]
        except ValueError:
            print "Unexpected format in xyz modification file:",self.modFile_,",line:",i+1
            sys.exit(2)
        # adjust indices if 1-indexed
        if self.indexed_ == 1:
            i0-=1
            i1-=1
        for xyzName in self.xyzNames_:
            myMol = pybel.readfile('xyz', xyzName).next()
            #First calculate the norm of the plane HO-O-H13
            print "Calculate the plane norm vector nv formed by HO-O-H11"
            coords = [ iatom.coords for iatom in myMol.atoms ]
            vOH= np.array(coords[iH])-np.array(coords[iO])
            vOH11= np.array(coords[iH11])-np.array(coords[iO])
            nv = np.cross(vOH,vOH11)
            nv = nv/np.linalg.norm(nv)
            #Then get the center point of the C1-C2 bond
            center = np.array(coords[iH11]) + vOH11/np.linalg.norm(vOH11)*bondLength*np.sqrt(3)/2.0
            #Get the C1-C2 direction
            c1c2 = np.cross(nv,vOH11)
            c1c2 = c1c2 / np.linalg.norm(c1c2)
            coordC1 = center + c1c2*bondLength/2.0
            coordC2 = center - c1c2*bondLength/2.0
            #Now translate C1,H12, H13 to the  new position of C1
            CHBondlength =1.1
            angle = np.radians(109.0)
            #v = coordC1 - coords[iC1]
            #coordH12 = coords[iH12] + v
            #coordH13 = coords[iH13] + v
            h =  coordC1 - coords[iH11] + coordC1-coordC2
            h = h/np.linalg.norm(h)
            centerH12H13 = coordC1 + h*CHBondlength*np.cos(angle/2.0)
            coordH12 = centerH12H13 + nv* CHBondlength*np.sin(angle/2.0)
            coordH13 = centerH12H13 - nv* CHBondlength*np.sin(angle/2.0)
            print "Attempt to move C(",iC1,")"
            atomToMove = myMol.OBMol.GetAtom(iC1+1)
            atomToMove.SetVector(pybel.ob.vector3(coordC1[0],coordC1[1],coordC1[2]))
            print "Attempt to move H(",iH12,")"
            atomToMove = myMol.OBMol.GetAtom(iH12+1)
            atomToMove.SetVector(pybel.ob.vector3(coordH12[0],coordH12[1],coordH12[2]))
            print "Attempt to move H(",iH13,")"
            atomToMove = myMol.OBMol.GetAtom(iH13+1)
            atomToMove.SetVector(pybel.ob.vector3(coordH13[0],coordH13[1],coordH13[2]))
            #Now translate C2,H21, H22 to the  new position of C2
            #v = coordC2 - coords[iC2]
            #coordH21 = coords[iH21] + v
            #coordH22 = coords[iH22] + v
            h =  coordC2 - coords[iH11] + coordC2-coordC1
            h = h/np.linalg.norm(h)
            centerH21H22 = coordC2 + h*CHBondlength*np.cos(angle/2.0)
            coordH21 = centerH21H22 + nv* CHBondlength*np.sin(angle/2.0)
            coordH22 = centerH21H22 - nv* CHBondlength*np.sin(angle/2.0)
            print "Attempt to move C(",iC2,")"
            atomToMove = myMol.OBMol.GetAtom(iC2+1)
            atomToMove.SetVector(pybel.ob.vector3(coordC2[0],coordC2[1],coordC2[2]))
            print "Attempt to move H(",iH21,")"
            atomToMove = myMol.OBMol.GetAtom(iH21+1)
            atomToMove.SetVector(pybel.ob.vector3(coordH21[0],coordH21[1],coordH21[2]))
            print "Attempt to move H(",iH22,")"
            atomToMove = myMol.OBMol.GetAtom(iH22+1)
            atomToMove.SetVector(pybel.ob.vector3(coordH22[0],coordH22[1],coordH22[2]))
            if moveAway == True:
                print "Finally move the ethylene away along the OH11 direction by 10 angstrom"
                v = vOH11/np.linalg.norm(vOH11)*10
                coordC1 = coords[iC1] + v
                coordH12 = coords[iH12] + v
                coordH13 = coords[iH13] + v
                coordC2 = coords[iC2] + v
                coordH21 = coords[iH21] + v
                coordH22 = coords[iH22] + v
                print "Attempt to move C(",iC1,")"
                atomToMove = myMol.OBMol.GetAtom(iC1+1)
                atomToMove.SetVector(pybel.ob.vector3(coordC1[0],coordC1[1],coordC1[2]))
                print "Attempt to move H(",iH12,")"
                atomToMove = myMol.OBMol.GetAtom(iH12+1)
                atomToMove.SetVector(pybel.ob.vector3(coordH12[0],coordH12[1],coordH12[2]))
                print "Attempt to move H(",iH13,")"
                atomToMove = myMol.OBMol.GetAtom(iH13+1)
                atomToMove.SetVector(pybel.ob.vector3(coordH13[0],coordH13[1],coordH13[2]))
                print "Attempt to move C(",iC2,")"
                atomToMove = myMol.OBMol.GetAtom(iC2+1)
                atomToMove.SetVector(pybel.ob.vector3(coordC2[0],coordC2[1],coordC2[2]))
                print "Attempt to move H(",iH21,")"
                atomToMove = myMol.OBMol.GetAtom(iH21+1)
                atomToMove.SetVector(pybel.ob.vector3(coordH21[0],coordH21[1],coordH21[2]))
                print "Attempt to move H(",iH22,")"
                atomToMove = myMol.OBMol.GetAtom(iH22+1)
                atomToMove.SetVector(pybel.ob.vector3(coordH22[0],coordH22[1],coordH22[2]))
            # delete the original xyz file
            subprocess.call(['rm','-r',xyzName])
            print "Successfully dissociated the C2H4",xyzName
            output=pybel.Outputfile("xyz",xyzName)
            output.write(myMol)
            output.close()
        return status

#TODO: Given a structure, add C2H4 group to the structure
#    def addC2H4(self):

    def modifyXyz(self):
        inp = open(self.modFile_).readlines()
        [mode,index]=["",0]
        try:
            [mode,index]= inp[0].strip('\n').split()
        except ValueError:
            print "Wrong format in line 1 of xyz modification file:", self.modFile_
        self.indexed_ = int(index)
        status = False
        if mode =="0":
             status=True
             #Don't do any modification.
        elif mode == "1":
           status=self.addHydrogens()
        elif mode == "2":
           status=self.moveHydrogens()
        elif mode == "3":
           status=self.modifyBonds()
        elif mode == "4": #Form the transition state of the reactions w/ loss of ethylene
           status=self.dissociateC2H4(moveAway=False)
        elif mode == "5":#Form the final state of the reactions w/ loss of ethylene
           status=self.dissociateC2H4(moveAway=True)
        elif mode == "6":#Delete the specified atoms
           status=self.removeAtoms()
        elif mode == "7":#Move a group of atoms along the direction given
           status=self.moveGroup()
        elif mode == "8":#rotate a set of atoms (set A) around the center of a plan defined by atoms (set B)
           status=self.rotateGroup()
        #elif mode == "9":#Add C2H42H4
        #   status=self.rotateGroup()
        else:
           print "Unhandled type of xyz modification"
           sys.exit(2)
        return status

    def generateModifiedXyz(self):
        print self.step_,'.Read in Xyz files and modify structures...'
        self.readConformers()
        status=self.modifyXyz()
        printLine()
        self.step_ +=1
        return True

   
#########################################
# 3) Generate input decks for TeraChem  #
#    to minimize the structures         #
#########################################
    def generateTcInput(self):
        print self.step_,'. Generating input decks...'
        #After doing all modifications, update the number of atoms /molecule
        self.numAtoms_ = pybel.readfile('xyz',self.xyzNames_[0]).next().OBMol.NumAtoms()
        for name,xyzname in itertools.izip(self.baseNames_,self.xyzNames_):
            #create working space for each conformer
            proc=subprocess.call(['mkdir','-p',name])
            #create input file
            proc=subprocess.call(['cp',xyzname, name])
            proc=subprocess.call(['cp','min.in',name])
            inputdeck=name+'/min.in'
            myfile=open(inputdeck,'a')
            myfile.write('coordinates '+name+'.xyz\n')
            myfile.write('end\n')
            myfile.close()
            #create job scheduling file
            script=name+'/sub.sh'
            data=open('submin.sh','r').read()
            data=data.replace('JOB_NAME',name)
            data=data.replace('XYZLENGTH',str(2+self.numAtoms_))
            target=open(script,'w')
            target.write(data)
            target.close()
        print 'Finished generating input decks!'
        printLine()
        self.step_ +=1
        
#######################################
# 4) Submit the minimization jobs and #
#    them to be finished. Keep track  #
#    conformer name vs energy         #
#######################################
    def checkJobs(self):
        # First Sleep and wait the squeue command to be able to display jobs 
        sleep(self.checkStatusTime_)
        while not jobFinished(self.jobIDs_, self.userName_):
            sleep(self.checkStatusTime_)
        print 'Finished minimizing all conformers!'
        
        #Now collect results
        print 'Checking minimization results...'
        subprocess.call('mkdir -p minimized',shell=True)
        for name in self.baseNames_:
            myOutputFile=name+'/min.out'
            myMinFile=name+'/min-'+name+'.xyz'
            cmd='grep Converged '+myOutputFile+' |wc -l'
            converged=int(subprocess.Popen(cmd,
                          shell=True,
                          stdout=subprocess.PIPE).communicate()[0].strip('\n'))
            if converged==0 :
               self.notConvergedNames_.append(name)
            else:
               subprocess.call(['cp',myMinFile,'minimized/'])
        
        if len(self.notConvergedNames_)>0:
           print 'Minimization of the following conformers are not converged:'
           for name in self.notConvergedNames_:
               print name
        else:
           print 'Minimization SUCCEEDED for all the conformers!'

    def submitJobs(self):
        print self.step_,'. Submit all the minimization jobs...'
        for name in self.baseNames_:
            subdir=self.workDir_+'/'+name
            os.chdir(subdir)
            proc=subprocess.Popen(['sbatch','sub.sh'],
                                  stdout=subprocess.PIPE)
            myJobID=proc.communicate()[0].strip('\n').split()[-1]
            self.jobIDs_.append(myJobID)
        
        os.chdir(self.workDir_)
        print 'Finished submitting jobs!\nWait for the jobs to be finished...\n'
        print 'My Slurm Job IDs:\n'
        for myJobID in self.jobIDs_:
            print myJobID

        self.checkJobs()
        printLine()
        self.step_ +=1
        
#######################################
# 5) Generate the folder of sorted    #
#    conformers                       #
#######################################
    def sortConformers(self):
        print self.step_,'. Sorting conformers by energy:'
        minFileList=subprocess.Popen(['ls','minimized'],stdout=subprocess.PIPE).communicate()[0].strip('\n').split()
        nameEnergyList=[]
        for xyzfile in minFileList:
           energy=readXyzEnergy('minimized/'+xyzfile)
           nameEnergyList.append([xyzfile,energy])
        
        nameEnergyList.sort(key=lambda tup:tup[1])
        #If the min_sorted directory already exists (for restart run), delete it
        if os.path.exists('min_sorted'):
           subprocess.call('rm -r min_sorted',shell=True)
        subprocess.call('mkdir -p min_sorted',shell=True)
        for i in range(0,len(nameEnergyList)):
            line=nameEnergyList[i]
            sourcename='minimized/'+line[0]
            targetname='min_sorted/'+'conf'+str(i)+'.xyz'
            print line[1],sourcename,'===>',targetname
            subprocess.call(['cp',sourcename,targetname])
        
        print 'Sorted conformers are stored in ./min_sorted'
        printLine()
        self.step_ +=1

################
# 6) Clean up  #  
################
    def writeCheckpoint(self):
        if os.path.exists('autoneb.chk'):
           subprocess.call('rm autoneb.chk',shell=True)
        if len(self.notConvergedNames_)>0:
           print 'Writing check point files for unconverged conformers.'
           output=open('autoneb.chk','w') 
           pickle.dump(self.notConvergedNames_, output)
           output.close()
            
# Delete the terachem run folders unless the conformer is not converged
    def cleanUp(self):
        '''TODO: check what needs to be modified for the restart mode'''
        print self.step_, '.Cleaning scratch files...'
        #if self.restart_==False:
        #   subprocess.call(['rm','-r',self.xyzDir_])
        
        for name in self.baseNames_:
            if not (name in self.notConvergedNames_):
               subprocess.call(['rm','-r',name])
        if len(self.notConvergedNames_) ==0 and self.no_sorting_ ==False:
           subprocess.call(['rm','-r','minimized'])
        self.writeCheckpoint()
        print 'AUTO NEB run finished!'
        printLine()

###############
#   restart   #
###############
    def resubmitJobs(self):
        print self.step_,'. Resubmit the unconverged minimization jobs...'
        self.baseNames_ = self.notConvergedNames_
        self.notConvergedNames_ = []
        if len(self.baseNames_) == 0:
            print 'No unconverged conformer exist. Nothing to do. Exiting...'
            return
        for name in self.baseNames_:
            subdir=self.workDir_+'/'+name
            os.chdir(subdir)
            # rename original structure
            xyzName = name+'.xyz'
            xyzlength = len(open(xyzName).readlines())
            # get last frame of optimized structure as the new starting point
            if os.path.exists('scr/optim.xyz'):
              subprocess.call(['mv',xyzName, 'original-'+xyzName])
              cmd = "tail -n "+ str(xyzlength)+" scr/optim.xyz > "+xyzName
              subprocess.call(cmd,shell=True)
            # submit unconverged jobs
            proc=subprocess.Popen(['sbatch','sub.sh'],
                                  stdout=subprocess.PIPE)
            myJobID=proc.communicate()[0].strip('\n').split()[-1]
            self.jobIDs_.append(myJobID)
        
        os.chdir(self.workDir_)
        print 'Finished submitting jobs!\nWait for the jobs to be finished...\n'
        print 'My Slurm Job IDs:\n'
        for myJobID in self.jobIDs_:
            print myJobID

        self.checkJobs()
        printLine()
        self.step_ +=1
         

#########################
# Driver for normal run #
#########################
    def runFromScratch(self):
        self.generateXyz()
        if self.no_opt_==False:
           self.generateTcInput()
           self.submitJobs()
           if self.no_sorting_ == False:
             self.sortConformers()
           self.cleanUp()

########################################
# Driver for run with xyz modification #
########################################
    def runWithModification(self):
        status=self.generateModifiedXyz()
        if status == True and self.no_opt_==False:
          self.generateTcInput()
          self.submitJobs()
          if self.no_sorting_ == False:
            self.sortConformers()
          self.cleanUp()
#########################
#   Driver for restart  #
#########################
    def restart(self): 
        self.resubmitJobs()
        if self.no_sorting_ == False:
          self.sortConformers()
        self.cleanUp()

#########################
#      Main Driver      #
#########################
    def run(self):
        self.readOpt()
        self.sanityCheck()
        if self.restart_ == True:
           self.restart()
        elif self.modFile_ != "":
           self.runWithModification()
        else:
           self.runFromScratch()
         
if __name__ == '__main__':
    gen = confgen()
    gen.run()
