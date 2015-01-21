README.txt

This project contains all the code required to recreate the results in the
Harvey 2015 paper on a parallel implementation of overlapping spheres
simulations: 'A parallel implementation of an off-lattice individual-based 
model of multicellular populations' by D. Harvey, A.G. Fletcher, J.M. Osborne 
and J.M. Pitt-Francis.


THIS PROJECT WILL ONLY WORK WITH CHASTE v3.2 (or v3.3).  The project was
developed with v3.2 but is known to work with v3.3.  It must be checked out
into a copy of the Chaste v3.2 (or v3.3) source code in

<chaste directory>/projects/Harvey2014

in order for file paths to be picked up correctly.

To run the examples you should do (for example):
cd <chaste directory>
scons build=GccOptNative projects/Harvey2014/test/TestUnitValidationLiteratePaper.hpp

There is a Chaste wiki page associated with this project that gives links to
more information on Chaste installation,
and the example simulations that are in the 'test' folder:
https://chaste.cs.ox.ac.uk/cgi-bin/trac.cgi/wiki/PaperTutorials/Harvey2014
