Python3 tools for to help preparing inlet velocity files for HemeLB simulations of branched domains.
Execution and brief explanation:

python3 InletBCPrep.py [inlets_radius.txt from voxelizer]
- returns inletsPU.txt which consists of inlet coordinates in PU and associated areas. This can be viewed in Paraview with the domain stl file to identify the location of iolets and indices

python3 VelocityFileGeneration.py [inlets_radius.txt from voxelizer]
- writes time-series files for each inlet [currently constant velocity, TODO scale a given heartbeat profile]
- specify the location and velocity of the true inlet, any known/desired flow rates at outlets controlled by velocity inlet files, and total flow rate going to pressure controlled outlets
- Flow rates at unspecified locations distributed either through equal flow at all sites or by area as a fraction total area of unknowns.

python3 weightsReader.py [filename.txt.weights.txt]
- reads weights file and returns statistics on number of points and average values
