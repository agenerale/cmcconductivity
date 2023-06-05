
import sys
from odbAccess import *
from abaqusConstants import *
from types import IntType
import numpy as np
import math
####################### READ IN FROM COMMAND LINE ##############################
# usage 'abaqus cae noGUI=avg_hflux.py -- '0001_JobQ'

num = sys.argv[-1]
temp = 1550
#odbName = num + '_' + str(temp[i]) +'_JobQ'
odbName = num + '_' + str(temp) +'_JobQ'
#print(odbName)
odb = openOdb(odbName + '.odb', readOnly=False)

################################ THICKNESS #################################

fileVoxel = num + '_Job.txt'
with open(fileVoxel, "rb") as f:
    lines = f.readlines()

linesnp = np.zeros(20) 
count = 0
for i in lines:
    if (count>0):
        linesnp[count-1]=float(lines[count])
    count += 1
        
height_max = linesnp[0]
height_min = linesnp[1]
width_max = linesnp[2]
width_min = linesnp[3]
xVoxel = int(linesnp[4])
yVoxel = int(linesnp[4])
zVoxel = int(linesnp[5])

tx = width_max - width_min
ty = width_max - width_min
tz = height_max - height_min

Tcold = 1500
Thot = 1600

################################ ELEMENT VOLUME #################################
steptoread=odb.steps['ThermalX']
frametoread=steptoread.frames[1]
odbEvol=frametoread.fieldOutputs['EVOL']
evolFieldValues=odbEvol.values
evol_vec = np.zeros(len(evolFieldValues)) 
evol_tot = 0
i = 0
for evol in evolFieldValues:
    evol_vec[i] = evol.data
    evol_tot += evol.data
    i += 1
  
allelements = odb.rootAssembly.instances['PART-1-1'].elementSets['ALLELEMENTS']   

################################## X DIRECTION ##################################
steptoread=odb.steps['ThermalX']
frametoread=steptoread.frames[1]
odbHflux=frametoread.fieldOutputs['HFL']
hflFieldValues=odbHflux.getSubset(region=allelements, position=CENTROID).values

hfl_vec = np.zeros(len(hflFieldValues))
i = 0
for hfl in hflFieldValues:
    hfl_vec[i] = ((hfl.data[0]**2)+(hfl.data[1]**2)+(hfl.data[2]**2))**0.5
    i += 1

hfl_x = 0
hfl_evol = np.zeros(len(evolFieldValues)) 
for i in range(evol_vec.shape[0]):
    hfl_evol[i] = hfl_vec[i]*(evol_vec[i]/evol_tot)
    hfl_x += hfl_evol[i]  

print 'hfl_x = %6.8f' % (hfl_x)

k11 = (hfl_x*tx)/(Thot-Tcold)
print 'k11 = %6.8f' % (k11)

################################## Y DIRECTION ##################################
steptoread=odb.steps['ThermalY']
frametoread=steptoread.frames[1]
odbHflux=frametoread.fieldOutputs['HFL']
hflFieldValues=odbHflux.getSubset(region=allelements, position=CENTROID).values

hfl_vec = np.zeros(len(hflFieldValues))
i = 0
for hfl in hflFieldValues:
    hfl_vec[i] = ((hfl.data[0]**2)+(hfl.data[1]**2)+(hfl.data[2]**2))**0.5
    i += 1

hfl_y = 0
hfl_evol = np.zeros(len(evolFieldValues)) 
for i in range(evol_vec.shape[0]):
    hfl_evol[i] = hfl_vec[i]*(evol_vec[i]/evol_tot)
    hfl_y += hfl_evol[i]   

print 'hfl_y = %6.8f' % (hfl_y)

k22 = (hfl_y*ty)/(Thot-Tcold)
print 'k22 = %6.8f' % (k22)

################################## Z DIRECTION ##################################
steptoread=odb.steps['ThermalZ']
frametoread=steptoread.frames[1]
odbHflux=frametoread.fieldOutputs['HFL']
hflFieldValues=odbHflux.getSubset(region=allelements, position=CENTROID).values

hfl_vec = np.zeros(len(hflFieldValues))
i = 0
for hfl in hflFieldValues:
    hfl_vec[i] = ((hfl.data[0]**2)+(hfl.data[1]**2)+(hfl.data[2]**2))**0.5
    i += 1

hfl_z = 0
hfl_evol = np.zeros(len(evolFieldValues)) 
for i in range(evol_vec.shape[0]):
    hfl_evol[i] = hfl_vec[i]*(evol_vec[i]/evol_tot)
    hfl_z += hfl_evol[i]   

print 'hfl_z = %6.8f' % (hfl_z)

k33 = (hfl_z*tz)/(Thot-Tcold)
print 'k33 = %6.8f' % (k33)

################################## WRITE OUT ##################################
data = file(odbName + '_k.txt', 'w')   
data.write('Thermal conductivity values in W/mK (K11, K22, K33)\n')
data.write(str(k11) + '\n')
data.write(str(k22) + '\n')
data.write(str(k33) + '\n')
data.close()  

odb.close()