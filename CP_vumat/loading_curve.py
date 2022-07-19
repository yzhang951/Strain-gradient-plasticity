# Yin Zhang , 17th Aug 2016
from odbAccess import *
import numpy as np
import math as mt
#import matplotlib as plt

#Reading output database

odbpath = 'Job-1.odb'

odb = openOdb(path=odbpath)

assembly = odb.rootAssembly

print 'Model data for ODB: ', odbpath

#Reading history output data

step = odb.steps['Step-1']

#print 'Node sets = ',odb.rootAssembly.nodeSets.keys()

#x_surface = assembly.nodeSets['_PickedSet18']

Node_Start = 1

Node_End = 20

# U = displacement, F = loading force

U = np.zeros(201)

F = np.zeros(201)

for i in range(Node_Start,Node_End+1,1):
    node = 'Node PART-1-1.' + str(i)
    region = step.historyRegions[node]
    u1data = region.historyOutputs['U1'].data
    f1data = region.historyOutputs['RF1'].data
    j = 0
    for time,force in f1data:
        F[j] = F[j] + force
        j = j + 1

j = 0
for time,disp in u1data:
    U[j] = disp      
    j = j + 1

dispFile = open('fitting.dat','w')

for i in range(201):
    j = i - 1; 
    dispFile.write("%10.4E  %10.4E\n" % ((U[j]/9)*100,F[j]/9))
dispFile.close()
