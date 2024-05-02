import numpy as np
import sys

path = sys.argv[1]

coordsList = []
areaList = []

dx=0.0001
shifts=np.array([-1593.39,-2302.05,-11046.6])

def transform_to_physical(pos, dx, shifts):
    return dx*(pos - shifts)

def between(a,low,high):
    return (a>low and a<high)

inletsPU=''

#path = "MESH0_INLET0_VELOCITY.txt.weights.txt"
print("starting to process", path)
f = open(path, "r")
for i in f:
    i = np.float_(i.split(','))

    coordsList.append([i[1],i[2],i[3]])
    areaList.append(i[5])
    c = transform_to_physical(np.array(coordsList[-1]),dx,shifts)
    inletsPU+= str(c[0])+" "+str(c[1])+" "+str(c[2])+" "+str(i[5])+"\n"
f.close()

with open("inletsPU.txt", "w") as outxml:
    outxml.write(inletsPU)
outxml.close()
