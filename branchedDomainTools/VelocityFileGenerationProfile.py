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
maxVelIN=0.02

knownRates=[(6,0.06),(8,0.05),(11,0.01),(21,0.01),(23,0.01),(12,0.03),(16,0.03),(59,0.01),(73,0.01),(69,0.01),(83,0.01),(78,0.01),(70,0.01),(68,0.01),(7,0.01),(4,0.01)]

outletRate=0.2

wavespeed=8.0 #m/s PWV is linked to cardiovascular disease and can vary widely depending on where measured and age

#Base profile, Flow_heartbeat.txt, has peak of 1m/s and cycle time of 1s (=60bpm)
#This script will generate a MESH0_INLET0_VELOCITY.txt file scaled to
# desired peakVelocity and period based on heartRate. If the simulation
# time extends beyond one cycle, it will cyclically repeat
#(from the 'back' for some reason, script accommodates this).

heartRate = 75 # rate in bpm
totalHeartBeats = 1 # total number of full heartbeats being simulated
#filename='InletProfile'+"_hRate"+str(heartRate)+"_MaxV"+str(peakVelocity)+".txt"
filename = 'INLET0_VELOCITY.txt'
filenameshift = 'INLET0_SHIFTVELOCITY.txt'
filenamedelay = 'INLET0_DELAYVELOCITY.txt'
warmupTime = 0.8 # time in sec of initialisation time




############
# Profile writing methods

def movingaverage(data, window_size):
    window = np.ones(int(window_size))/float(window_size)
    #return np.convolve(data, window, 'valid')
    return np.concatenate((data[0], np.convolve(data, window, 'valid'), data[-1]), axis=None)

def movave2(data, window_size):
    cumsum_vec = numpy.cumsum(numpy.insert(data, 0, 0))
    return (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width

def refine(data,pts):
    X2=[]
    for i in range(len(data)-1):
        X2 = X2 + list(np.linspace(data[i],data[i+1],pts,endpoint=False))
    return X2 + [data[-1]]

def generateTotalProfile(baseProfile):
    f = open(baseProfile, "r")
    X=[]
    Y=[]
    Xa = []
    Ya = []
    Xr = []
    Yr = []

    for i in f:
        i = np.float_(i.split())
        if np.shape(i)==(0,):
            break
        X.append(i[0])
        Y.append(i[1])

    f.close()

    Xr = X
    Yr = Y
    for i in range(2):
        Xr = refine(movingaverage(Xr,8),10)
        Yr = refine(movingaverage(Yr,8),10)

    maxR = max(Xr)
    Xr = [60*Xr[i]/(heartRate*maxR) for i in range(len(Xr))] #scale for rate
    maxR = max(Yr)
    Yr = [Yr[i]/maxR for i in range(len(Yr))] #scale for Velocity

    Yf = Yr

    Xa.append(0.0); Ya.append(0.0)
    Xa.append(0.25*warmupTime); Ya.append(Yf[0])
    for j in range(totalHeartBeats):
        for i in range(len(Xr)):
            Xa.append(warmupTime + Xr[i] + j*max(Xr))
            Ya.append(Yf[i])

    return Xa, Ya

def writeScaledProfile(filename,maxVel,area,baseX,baseY):
    print("mean velocity ", 0.5*maxVel*np.average(baseY), "m/s")
    print("mean flow rate ", 0.5*maxVel*np.average(baseY)*area*6e7, "mL/min")

    maxR = max(baseY)
    Yr = [maxVel*baseY[i]/maxR for i in range(len(baseY))] #scale for Velocity

    contents='0.0 0.0\n'
   
    for i in range(1,len(baseX)):
        contents += str(baseX[i]) + ' ' + str(Yr[i]) + '\n'
    
    f = open(filename,'w')
    f.write(contents)
    f.close()

    return 

def rotate(li, x):
  return li[-x % len(li):] + li[:-x % len(li)]

#Delays profile by extending the initialisation period
def writeScaledDelayedProfile(filename,maxVel,area,delay,baseX,baseY):
    print("mean velocity ", 0.5*maxVel*np.average(baseY), "m/s")
    print("mean flow rate ", 0.5*maxVel*np.average(baseY)*area*6e7, "mL/min")

    maxR = max(baseY)
    Yr = [maxVel*baseY[i]/maxR for i in range(len(baseY))] #scale for Velocity
   
    print('delay of', delay, 'sec')
    contents='0.0 0.0\n'
    contents += str(baseX[1]) + ' ' + str(Yr[1]) + '\n'
   
    for i in range(2,len(baseX)):
        contents += str(baseX[i]+delay) + ' ' + str(Yr[i]) + '\n'
    
    f = open(filename,'w')
    f.write(contents)
    f.close()

    return 

#Delays profile by shifting the whole heartbeat portion by delay period
def writeScaledShiftedProfile(filename,maxVel,area,delay,baseX,baseY):
    print("mean velocity ", 0.5*maxVel*np.average(baseY), "m/s")
    print("mean flow rate ", 0.5*maxVel*np.average(baseY)*area*6e7, "mL/min")

    maxR = max(baseY)
    Yr = [maxVel*baseY[i]/maxR for i in range(len(baseY))] #scale for Velocity
    print('delay of', delay, 'sec')
    print('shift by ', np.floor((60*delay/heartRate)*len(Yr[2:])/totalHeartBeats).astype(int), ' steps')

    Yb = rotate(Yr[2:], np.floor((60*delay/heartRate)*len(Yr[2:])/totalHeartBeats).astype(int))

   
    contents='0.0 0.0\n'
    contents += str(baseX[1]) + ' ' + str(Yb[0]) + '\n'
   
    for i in range(2,len(baseX)):
        contents += str(baseX[i]) + ' ' + str(Yb[i-2]) + '\n'
    
    f = open(filename,'w')
    f.write(contents)
    f.close()

    return 

#######################
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

qIN = 0.5*maxVelIN*areaList[trueInletIndex]
areaSum = sum(areaList) - areaList[trueInletIndex] - sum(areaList[j] for j in knownLets)
predictedFractions = ''

bX, bY = generateTotalProfile("Flow_heartbeat.txt")

for i in range(len(areaList)):
    delay=dx*np.linalg.norm(np.array(coordsList[i])-np.array(coordsList[trueInletIndex]))/wavespeed 
    if i==trueInletIndex:
        writeScaledProfile(filename.replace('0',str(i)),maxVelIN,areaList[i],bX,bY)
        writeScaledDelayedProfile(filenamedelay.replace('0',str(i)),maxVelIN,areaList[i],delay,bX,bY)
        writeScaledShiftedProfile(filenameshift.replace('0',str(i)),maxVelIN,areaList[i],delay,bX,bY)
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(1.0) + ' ' + str(areaList[i]) + '\n'
    elif i in knownLets:
        vmax = 2.0*knownRates[knownLets.index(i)][1]*qIN/areaList[i]
        
        if vmax>globalVMax:
            globalVMax = vmax

        print(i)
        writeScaledProfile(filename.replace('0',str(i)),-1.0*vmax,areaList[i],bX,bY)
        writeScaledDelayedProfile(filenamedelay.replace('0',str(i)),-1.0*vmax,areaList[i],delay,bX,bY)
        writeScaledShiftedProfile(filenameshift.replace('0',str(i)),-1.0*vmax,areaList[i],delay,bX,bY)
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(knownRates[knownLets.index(i)][1]) + ' ' + str(areaList[i]) + '\n'
    else:
        #vmax = 2.0*(1.0 - outletRate - totalKnown)*qIN/(areaList[i]*(len(areaList)-len(knownLets)-1)) #even distribution of flow to unassigned outlets
        vmax = 2.0*(1.0 - outletRate - totalKnown)*qIN/(areaSum) #distribution of flow to unassigned outlets based on area

        if vmax>globalVMax:
            globalVMax = vmax
        print(i)
        writeScaledProfile(filename.replace('0',str(i)),-1.0*vmax,areaList[i],bX,bY) 
        writeScaledDelayedProfile(filenamedelay.replace('0',str(i)),-1.0*vmax,areaList[i],delay,bX,bY) 
        writeScaledShiftedProfile(filenameshift.replace('0',str(i)),-1.0*vmax,areaList[i],delay,bX,bY) 
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(0.5*vmax*areaList[i]/qIN) + ' ' + str(areaList[i]) + '\n'

with open("PredictedFlowFractions.txt", "w") as outxml:
    outxml.write(predictedFractions)
outxml.close()
print("Global VMax = ", globalVMax, ", for estimated Ma of ", 3*globalVMax*dt/dx)
