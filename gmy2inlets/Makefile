# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE. 

# For Copyright and Licensing please refer to the LICENSE and
# THIRD_PARTY_LICENSES file in the top level directory of this package

#In general, usr may refer to HemePure/dep/install directory
#CFLAGS = -I/<folderlocation>/HemePure/dep/install/include/tirpc
#LFLAGS = -L/<folderlocation>/HemePure/dep/install/lib

all: gmy2lets sgmy2lets

debug: gmy2lets_debug sgmy2lets_debug

clean:
	rm -f gmy2lets sgmy2lets
	
gmy2lets: gmy2lets.cc
	g++ -O3 -Wall -I/usr/include/tirpc gmy2lets.cc -lz -L/usr/lib -ltirpc -o gmy2lets

gmy2lets_debug: gmy2lets.cc
	g++ -Wall -g -I/usr/include/tirpc gmy2lets.cc -lz -L/usr/lib -ltirpc -o gmy2lets

sgmy2lets: sgmy2lets.cc
	g++ -O3 -Wall -I/usr/include/tirpc sgmy2lets.cc -lz -L/usr/lib -ltirpc -o sgmy2lets

sgmy2lets_debug: sgmy2lets.cc
	g++ -Wall -g -I/usr/include/tirpc sgmy2lets.cc -lz -L/usr/lib -ltirpc -o sgmy2lets

