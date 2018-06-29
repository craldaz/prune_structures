import subprocess
import sys
import copy
import os
from time import sleep
import getopt
import pybel
import numpy as np
import math
import pickle
import rotate
import itertools
#input options
xyzName1=sys.argv[1]
xyzName2=sys.argv[2]

def getCoords(mymol):
    return np.array([ iatom.coords for iatom in mymol.atoms ]) 

def asignCoords(coords,myMol):
  for i in range(0,len(coords)):
       myMol.OBMol.GetAtom(i+1).SetVector(pybel.ob.vector3(coords[i][0],coords[i][1],coords[i][2]))

def getXYZData(myMol):
    coords=getCoords(myMol)
    anums=getAtomicNums(myMol)
    xyzData=[]
    for i in range(0,len(coords)):
       xyzData.append([anums[i],coords[i][0],coords[i][1],coords[i][2]])
    return xyzData

def getAtomicNums(mymol):
    return [ iatom.atomicnum for iatom in mymol.atoms ]

#Use the phi,theta,psi convention by Golstein and Landau (Walfram Language)
# http://mathworld.wolfram.com/EulerAngles.html
'''
   phi,psi: [-pi,pi]
   theta: [-pi/2,pi/2]
'''
def rotateEuler(phi,theta,psi,coords):
    R=np.zeros((3,3))
    rPhi=np.deg2rad(phi)
    rTheta=np.deg2rad(theta)
    rPsi = np.deg2rad(psi)
    cosPsi=np.cos(rPsi)
    sinPsi=np.sin(rPsi)
    cosTheta=np.cos(rTheta)
    sinTheta=np.sin(rTheta)
    cosPhi=np.cos(rPhi)
    sinPhi=np.sin(rPhi)
    R[0][0] = cosPsi*cosPhi - cosTheta*sinPhi*sinPsi
    R[0][1] = cosPsi*sinPhi + cosTheta*cosPhi*sinPsi
    R[0][2] = sinPsi*sinTheta
    R[1][0] = -sinPsi*cosPhi - cosTheta*sinPhi*cosPsi
    R[1][1] = -sinPsi*sinPhi + cosTheta*cosPhi*cosPsi
    R[1][2] = cosPsi*sinTheta
    R[2][0] = sinTheta*sinPhi
    R[2][1] = -sinTheta*cosPhi
    R[2][2] = cosTheta
    print R
    newcoords=np.transpose(np.dot(R,np.transpose(coords)))
    return newcoords

def saveXyz(myMol,xyzname):    
    output=pybel.Outputfile("xyz",xyzname)
    output.write(myMol)
    output.close()

def getRecMatrix(xyzData):
    n =len(xyzData) 
    M=np.zeros((n,n))
    for i in range(0,n):
       for j in range(i+1,n):
         dist = np.linalg.norm(np.array(xyzData[i][1:])-np.array(xyzData[j][1:]))
         if dist<1e-6:
             print 'divide by 0!:',i,j,dist
         M[i][j]=xyzData[i][0]*xyzData[j][0]/dist
         M[j][i]=M[i][j]
    return M
      

def getRMSDRecMatrix(xyzData1,xyzData2):
    coords1=np.array([ xyzData1[i][1:] for i in range(0,len(xyzData1)) ])
    coords2=np.array([ xyzData2[i][1:] for i in range(0,len(xyzData2)) ])
    rec1=getRecMatrix(xyzData1)
    rec2=getRecMatrix(xyzData2)
    dist=np.linalg.norm(rec1-rec2)
    n = rec1.shape[0]
    rmsd=dist/n
    return rmsd

def printXyz(xyzData):
    output=open('relabel-'+xyzName2,'w')
    output.write(str(len(xyzData))+'\n\n')
    for line in xyzData:
        if line[0] == 1:
            output.write('H ')
        elif line[0] == 6:
            output.write('C ')
        elif line[0] == 8:
            output.write('O ')
        elif line[0] == 7:
            output.write('N ')
        output.write(str(line[1])+' ')
        output.write(str(line[2])+' ')
        output.write(str(line[3])+'\n')
             

# SetDict is a dictionary that maps atomic numbers to sets containing indices of atoms with that atomic number
def getElementSets(anums):
    setDict=dict()
    for i in range(0,len(anums)):
         if anums[i] not in setDict:
              setDict[anums[i]]=set()
         setDict[anums[i]].add(i)
    return setDict

def permutateOneCoord(xyzData, newId, oldId):
    #print xyzData[oldId]
    #print xyzData[newId]
    newData=copy.deepcopy(xyzData)
    for i in range(0,4):
        newData[newId][i]=xyzData[oldId][i]
        newData[oldId][i]=xyzData[newId][i]
    #print xyzData[oldId]
    #print xyzData[newId]
    #print newData[newId]
    #print newData[oldId]
    #print 'permutating ',newId,' with ',oldId
    return newData

def permutateOneCoordInPlace(xyzData, newId, oldId):
    tmp = [ x for x in xyzData[newId]]
    for i in range(0,4):
        xyzData[newId][i]=xyzData[oldId][i]
    for i in range(0,4):
        xyzData[oldId][i]=tmp[i]


def autoRelabel(xyzData1, xyzData2):        
    print 'Looking for pertumation that make the labeling consistent'
    n=len(xyzData1)
    changed_counter=10
    while changed_counter > 0:
       changed_counter = 0
       for i in range(0,n):
           an1=xyzData1[i][0]
           distdict=dict()
           min_dist=100000
           min_id = -1
           for j in range(i,n):
               an2=xyzData2[j][0]
               if an2==an1:
                   newXyzData2=permutateOneCoord(xyzData2,i,j)     
                   distdict[j]=getRMSDRecMatrix(xyzData1,newXyzData2) 
           for j in distdict:
               if distdict[j] < min_dist:
                    min_dist=distdict[j]
                    min_id=j
           if i != min_id:
              changed_counter +=1
              permutateOneCoordInPlace(xyzData2,i,min_id)
       print changed_counter 
            
       
# Here start the main funcitons
myMol1 = pybel.readfile('xyz', xyzName1).next()
myMol2 = pybel.readfile('xyz', xyzName2).next()
xyzData1=getXYZData(myMol1)
xyzData2=getXYZData(myMol2)

autoRelabel(xyzData1,xyzData2)
printXyz(xyzData2)
         
