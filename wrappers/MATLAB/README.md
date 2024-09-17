# MATLAB bindings


These bindings are still WIP. `DAGGER` uses a wrapper class to convert MATLAB arrays to `std::vector`-like data structure (that works) and the other way (that does not work yet).

Relevant files:

- `./matmain.hpp` -> contains a minimal example using the main classes to be binded with `clibgen` (basically if it works, I can make everything work)
- `./builder.m` -> MATLAB script to build the minimal example (Adapted to linux, potential for windows too)
- `./demo/tester.m` -> runs a simple test once built

Tested on `Ubuntu 22.04` with `g++11`. Note that I need to run MATLAB in a shell script and edit the `LD_LIBRARY_PATH` to get the building to work:
```
export LD_LIBRARY_PATH=/home/bgailleton/MATLAB/extern/bin/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/bgailleton/MATLAB/extern/bin:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/bgailleton/MATLAB/bin/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/bgailleton/MATLAB/bin:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/bgailleton/Desktop/code/DAGGER/wrappers/MATLAB/DAGGER:$LD_LIBRARY_PATH
```
I also had to change the `libstdc++.so.6` of MATLAB folders to match my own `g++11`. On windows, it only worked with the `MinGW++` extention from MathWork.























<!--  -->
