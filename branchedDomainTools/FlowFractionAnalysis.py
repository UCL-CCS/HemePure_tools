import numpy as np
import sys

predicted = sys.argv[1]
observed = sys.argv[2]

coordsList = []
areaList = []
fractionList = []


def transform_to_physical(pos, dx, shifts):
    return dx*(pos - shifts)

def between(a,low,high):
    return (a>low and a<high)


####################
#Set these values to your particular requirements
dx=0.0001
dt=0.00008
shifts=np.array([-1593.39,-2302.05,-11046.6])

trueInletIndex=5
maxVelIN=0.01

knownRates=[(6,0.2),(7,0.06),(10,0.05),(65,0.01),(51,0.01),(60,0.01),(52,0.01),(50,0.01),(55,0.01),(42,0.01)]

outletRate=0.2
############

totalKnown = sum(j for i, j in knownRates)
knownLets=[i for i,j in knownRates]
globalVMax = maxVelIN


print("starting to process", predicted)
f = open(predicted, "r")
for i in f:
    i = np.float_(i.split(' '))

    coordsList.append([i[0],i[1],i[2]])
    fractionList.append(i[3])
    areaList.append(i[4])

f.close()

velSum = [0.0 for i in range(len(areaList))]
velCount = [0 for i in range(len(areaList))]

print("starting to process", observed)
f = open(observed, "r")
next(f)
for i in f:
    i = np.float_(i.split(' '))
    dist=1.0e10
    for j in range(len(areaList)):
        d = np.sqrt((i[1]-coordsList[j][0])**2 +(i[2]-coordsList[j][1])**2 +(i[3]-coordsList[j][2])**2)
        if d < dist:
            dist = d
            close = j;
    velSum[close] += np.sqrt(i[4]**2 + i[5]**2 + i[6]**2)
    velCount[close] += 1

f.close()

qIN = areaList[trueInletIndex]*velSum[trueInletIndex]/velCount[trueInletIndex]
for i in range(len(areaList)):
    #print('index ', i, ' ratio ', (areaList[i]*velSum[i]/velCount[i])/(0.5*maxVelIN*areaList[trueInletIndex] * fractionList[i]))
    #print('index ', i, ' ratio ', (areaList[i]*velSum[i]/velCount[i])/(qIN * fractionList[i]))
    print('inlet index ', i, ' fraction ', (areaList[i]*velSum[i]/velCount[i])/qIN)
print('total fraction', sum( (areaList[i]*velSum[i]/velCount[i])/qIN for i in range(len(areaList)))-1.0 )

#########################
#Outlets
#velSum = [0.0 for i in range(len(areaList))]
#velCount = [0 for i in range(len(areaList))]
coordsList=[]
areaList=[]

path='outlets_radius.txt'
print('starting to process ', path)
f = open(path, "r")
for i in f:
    i = np.float_(i.split(','))

    coordsList.append([dx*i[1],dx*i[2],dx*i[3]])
    areaList.append(i[5])
f.close()

velSum = [0.0 for i in range(len(areaList))]
velCount = [0 for i in range(len(areaList))]

observed=observed.replace('in','out')
print("starting to process", observed)
f = open(observed, "r")
next(f)
for i in f:
    i = np.float_(i.split(' '))
    dist=1.0e10
    for j in range(len(areaList)):
        d = np.sqrt((i[1]-coordsList[j][0])**2 +(i[2]-coordsList[j][1])**2 +(i[3]-coordsList[j][2])**2)
        if d < dist:
            dist = d
            close = j;
    velSum[close] += np.sqrt(i[4]**2 + i[5]**2 + i[6]**2)
    velCount[close] += 1

f.close()

for i in range(len(areaList)):
    #print('index ', i, ' ratio ', (areaList[i]*velSum[i]/velCount[i])/(0.5*maxVelIN*areaList[trueInletIndex] * fractionList[i]))
    #print('index ', i, ' ratio ', (areaList[i]*velSum[i]/velCount[i])/(areaList[trueInletIndex]*velSum[trueInletIndex]/velCount[trueInletIndex] * fractionList[i]))
    print('outlet index ', i, ' fraction ', (areaList[i]*velSum[i]/velCount[i])/qIN)
print('total fraction', sum( (areaList[i]*velSum[i]/velCount[i])/qIN for i in range(len(areaList))) )

