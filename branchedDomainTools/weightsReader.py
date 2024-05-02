import numpy as np
import sys

path = sys.argv[1]
#path = "MESH0_INLET0_VELOCITY.txt.weights.txt"
print("starting to process", path)
f = open(path, "r")

wtSum = 0
wtCount = 0
zeroCount = 0
wtMax = 0
wtMin = 10
wtMinNZ = 10

for i in f:  
    i = np.float_(i.split())

    wtSum = wtSum + i[3]
    wtCount = wtCount + 1
    
    if i[3] == 0.0:
        zeroCount = zeroCount + 1

    if ((i[3] < wtMinNZ) and i[3]>0.0):
        wtMinNZ = i[3]
    
    if i[3] < wtMin:
        wtMin = i[3]
        
    if i[3] > wtMax:
        wtMax = i[3]

print("For ", wtCount, " points:")
print(zeroCount, " are wt=0.0 (=", zeroCount/wtCount*100,"%)")
print("Max = ", wtMax, ", Min = ", wtMin, ", Min non-zero = ", wtMinNZ, " Ave = ", wtSum/wtCount, ", non-zero Ave = ", wtSum/(wtCount-zeroCount))

