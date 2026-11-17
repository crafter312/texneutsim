This README is for the texneutsim Geant4 simulation.

To use this code after downloading, starting from project root directory:
$ git submodule update --init --recursive
$ cmake -S . -B build
$ cd build
$ make -j

To clean build, starting from project root directory:
$ cd build
$ make clean
$ make -j

To super duper clean build, starting from project root directory:
$ rm -rf build
$ cmake -S . -B build
$ cd build
$ make -j

To manage git submodules, I prefer to do so manually from inside
each submodule. Each submodule is treated as its own git repository,
so you can just CD into, for example, `lib/Li6sim/` and manage it
with all the standard git commands (i.e. `git status`, `git pull`, even
`git push` if you make changes to the submodule in place).

Now that the package is built we can run it from the build directory. To run
there are several options. The first is to run it as is with the visualization.
To do this, simply do
$ ./texneutsim

This will execute the code with the vis.mac macro file. To run it with any 
specific macro file, use your intended macro file as an argument:
$ ./texneutsim run_6Li.mac
