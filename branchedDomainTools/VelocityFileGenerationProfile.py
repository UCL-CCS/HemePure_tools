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


#Base profile, Flow_heartbeat.txt, has peak of 1m/s and cycle time of 1s (=60bpm)
#This script will generate a MESH0_INLET0_VELOCITY.txt file scaled to
# desired peakVelocity and period based on heartRate. If the simulation
# time extends beyond one cycle, it will cyclically repeat
#(from the 'back' for some reason, script accommodates this).
peakVelocity = 0.4/1.4 #Ultimate maximum velocity [m/s] (throat factor for fistula arteries)

heartRate = 75 # rate in bpm
totalHeartBeats = 1 # total number of full heartbeats being simulated
#filename='InletProfile'+"_hRate"+str(heartRate)+"_MaxV"+str(peakVelocity)+".txt"
filename = 'INLET0_VELOCITY.txt'
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

def writeScaledProfile(filename,maxVel,area,baseProfile):
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


    print("mean velocity ", 0.5*maxVel*np.average(Yr), "m/s")
    print("mean flow rate ", 0.5*maxVel*np.average(Yr)*area*6e7, "mL/min")

    maxR = max(Xr)
    Xr = [60*Xr[i]/(heartRate*maxR) for i in range(len(Xr))] #scale for rate
    maxR = max(Yr)
    Yr = [maxVel*Yr[i]/maxR for i in range(len(Yr))] #scale for Velocity

    Yf = Yr

    f = open(filename,'w')
    contents = '0.0 0.0\n'
    contents += str(0.25*warmupTime) + ' ' + str(Yf[0]) + '\n'

    Xa.append(0.0); Ya.append(0.0)
    Xa.append(0.25*warmupTime); Ya.append(Yf[0])
    for j in range(totalHeartBeats):
        for i in range(len(Xr)):
            contents += str(warmupTime + Xr[i] + j*max(Xr)) + ' ' + str(Yf[i]) + '\n'
            Xa.append(warmupTime + Xr[i] + j*max(Xr))
            Ya.append(Yf[i])
    f.write(contents)
    f.close()

    #plt.plot(Xa,Ya,'.-')
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

for i in range(len(areaList)):
    
    if i==trueInletIndex:
        writeScaledProfile(filename.replace('0',str(i)),maxVelIN,areaList[i],"Flow_heartbeat.txt")
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(1.0) + ' ' + str(areaList[i]) + '\n'
    elif i in knownLets:
        vmax = 2.0*knownRates[knownLets.index(i)][1]*qIN/areaList[i]
        
        if vmax>globalVMax:
            globalVMax = vmax

        print(i,vmax )
        writeScaledProfile(filename.replace('0',str(i)),-1.0*vmax,areaList[i],"Flow_heartbeat.txt")
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(knownRates[knownLets.index(i)][1]) + ' ' + str(areaList[i]) + '\n'
    else:
        #vmax = 2.0*(1.0 - outletRate - totalKnown)*qIN/(areaList[i]*(len(areaList)-len(knownLets)-1)) #even distribution of flow to unassigned outlets
        vmax = 2.0*(1.0 - outletRate - totalKnown)*qIN/(areaSum) #distribution of flow to unassigned outlets based on area

        if vmax>globalVMax:
            globalVMax = vmax
        print(i,vmax)
        writeScaledProfile(filename.replace('0',str(i)),-1.0*vmax,areaList[i],"Flow_heartbeat.txt") 
        predictedFractions+=str(coordsList[i][0]*dx) + ' '+str(coordsList[i][1]*dx) + ' '+str(coordsList[i][2]*dx) + ' '+str(0.5*vmax*areaList[i]/qIN) + ' ' + str(areaList[i]) + '\n'

with open("PredictedFlowFractions.txt", "w") as outxml:
    outxml.write(predictedFractions)
outxml.close()
print("Global VMax = ", globalVMax, ", for estimated Ma of ", 3*globalVMax*dt/dx)
