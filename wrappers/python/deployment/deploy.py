'''
This scripts deploys the
'''

import os
import shutil
import subprocess as sub
import glob



# Find all files in the source directory with the ".txt" extension
files_to_copy = glob.glob("../../../DAGGER/*.hpp")
files_to_copy += glob.glob("../../../fastflood/*.hpp")
files_to_copy += glob.glob("../../../popscape/*.hpp")
files_to_copy += glob.glob("../*.hpp")

print("Copying all the files ... ")

for file in files_to_copy:
	shutil.copy(file, "./dagger/includes/")

shutil.copy("../main.cpp", "./dagger/main.cpp")
print("Done!")
