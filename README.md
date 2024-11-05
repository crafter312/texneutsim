This simulation is used to simulate the TexNeut neutron detector array with
scintillation light. 

To use this code, after downloading: 
$ cd build	
$ cmake ..
c to configure
make sure paths are appropriate for your install of Geant4
c to configure again
assuming no errors, g to generate
q to quit
$ make clean
$ make -jN

Now that the package is built we can run it from the build directory. To run
there are several options. The first is to run it as is with the visualization.
To do this,
simply do
$ ./examplesim
This will execute the code with the vis.mac macro file. To run it with any 
specific macro file, use your intended macro file as the second argument with
after the examplesim execution, i.e., 
$ ./examplesim run.mac

In this branch we will add radioactive decay and a scintillation type material