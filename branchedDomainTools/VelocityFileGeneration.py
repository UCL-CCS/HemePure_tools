import numpy as np
import sys

path = sys.argv[1]

coordsList = []
areaList = []


def transform_to_physical(pos, dx, shifts):
    return dx*(pos - shifts)

def between(a,low,high):
    return (a>low and a<high)

inletsPU=''

####################
#Set these values to your particular requirements
dx=0.0001
dt=0.00008

trueInletIndex=5
maxVelIN=0.01

knownRates=[(6,0.2),(7,0.06),(10,0.05),(65,0.01),(51,0.01),(60,0.01),(52,0.01),(50,0.01),(55,0.01),(42,0.01)]

outletRate=0.2
############

totalKnown = sum(j for i, j in knownRates)
knownLets=[i for i,j in knownRates]
globalVMax = maxVelIN

#path = "MESH0_INLET0_VELOCITY.txt.weights.txt"
print("starting to process", path)
f = open(path, "r")
for i in f:
    i = np.float_(i.split(','))

    coordsList.append([i[1],i[2],i[3]])
    areaList.append(i[5])

f.close()

qIN = maxVelIN*areaList[trueInletIndex]
areaSum = sum(areaList) - areaList[trueInletIndex] - sum(areaList[j] for j in knownLets)
predictedFractions = ''

for i in range(len(areaList)):
    profile="0.0 0.0\n"
    
    if i==trueInletIndex:
        profile+="0.5 "+str(maxVelIN)+"\n10.0 "+str(maxVelIN)
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(1.0) + '\n'
    elif i in knownLets:
        vmax = 2.0*knownRates[knownLets.index(i)][1]*qIN/areaList[i]
        
        if vmax>globalVMax:
            globalVMax = vmax

        print(i,vmax )
        profile+="0.5 "+str(-1.0*vmax)+"\n10.0 "+str(-1.0*vmax)
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(knownRates[knownLets.index(i)][1]) + '\n'
    else:
        #vmax = 2.0*(1.0 - outletRate - totalKnown)*qIN/(areaList[i]*(len(areaList)-len(knownLets)-1)) #even distribution of flow to unassigned outlets
        vmax = 2.0*(1.0 - outletRate - totalKnown)*qIN/(areaSum) #distribution of flow to unassigned outlets based on area

        if vmax>globalVMax:
            globalVMax = vmax
        print(i,vmax)
        profile+="0.5 "+str(-1.0*vmax)+"\n10.0 "+str(-1.0*vmax)
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(0.5*vmax*areaList[i]/qIN) + '\n'
    with open("INLET"+str(i)+"_VELOCITY.txt", "w") as outxml:
        outxml.write(profile)
    outxml.close()

with open("PredicedFlowFractions.txt", "w") as outxml:
    outxml.write(predictedFractions)
outxml.close()
print("Global VMax = ", globalVMax, ", for estimated Ma of ", 3*globalVMax*dt/dx)
