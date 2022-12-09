% Script building the DAGGER bindings for MATLAB
% WARNING: WIP, not really stable yet
% So far only works with MinGW64 installed with MATLAB WIndows
% will probably work on linux, I just need to fix a few thing
% As of the 6/12/2022 it does not work on MacOS: the official example to bind c++ and MATLAB does not even work...

% Locating hte header(s) of interest (file holding the direct matlab/c++ adaptation)
hppfiles = ["matmain.hpp"];

% Path required to load the MATLAB and DAGGER-related inputs
includepath = [fullfile(matlabroot,"extern","include"), "../../DAGGER"];

% Libraries enabling management of MATLAB array, vectors and matrices
library_files = [fullfile(matlabroot,"extern","lib","win64","mingw64","libMatlabDataArray.lib"),...
    fullfile(matlabroot,"extern","lib","win64","mingw64","libMatlabEngine.lib")];

% clibgen generates the actual binding files, all the magic happens here
% MATLAB nicely automates the thing, but getting to this point was painful
% and debugging messages are very, very unfriendly
clibgen.generateLibraryDefinition(hppfiles,... % Headers
    "OverwriteExistingDefinitionFiles",true, ... % Yes I want to reload
    "DefinedMacros","DAGGER_FT_MATLAB", ... % DEfines the compiling macro enabling data conversion for MATLAB in the DAGGER code (each language has its own)
    "Libraries",library_files, ... % Uses the MATLAB libraries (.lib linked to MinGW)
    IncludePath=includepath,... % finds the connected header files to include on the c++ side
    AdditionalCompilerFlags=["-std=c++17" "-O3"],... % Compiling options: I use standard c++17 constructs, and optimise to level 3
    PackageName="DAGGER"); % in clib, the package will be ingested as DAGGER

build(defineDAGGER);