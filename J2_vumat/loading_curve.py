# Yin Zhang , 17th Aug 2016
from odbAccess import *
import numpy as np
import math as mt
#import matplotlib as plt

#Reading output database

odbpath = 'J_2000.odb'

odb = openOdb(path=odbpath)

assembly = odb.rootAssembly

print 'Model data for ODB: ', odbpath 
#Reading history output data

step = odb.steps['Step-1']

print 'Node sets = ',odb.rootAssembly.nodeSets.keys()

#x_surface = assembly.nodeSets['_PickedSet18']

Node_Start = 6

Node_End = 216

# U = displacement, F = loading force

U = np.zeros(12000)

F = np.zeros(12000)

for i in range(1,162,1):
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

for i in range(300):
    j = i - 1; 
#   dispFile.write("%10.4E  %10.4E\n" % (mt.log1p(U[j])*100,(1+U[j])*F[j]/1E6))
    dispFile.write("%10.4E  %10.4E\n" % (U[j]*100/0.2,F[j]/1E6/16/0.2))
dispFile.close()
