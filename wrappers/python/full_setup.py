'''
This script runs a full install of the dagger-related codes
B.G.
'''

import subprocess as sub
import sys

print("\n\n~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=")
print("I am going to install/update the full dagger suite in python.\n\n")
# print("Note: it automates the full process to force the update and removes all the build/dist subfolders (if you don't know what it is just ignore)")

if("-v" in sys.argv):
	stttd = None
else:
	print("use option -v to get the whole install log")
	stttd = sub.DEVNULL

print("Installing dagger ... ", end = "", flush = True)
try:
	sub.run("rm -r dist", cwd = "./dagger/", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("rm -r build", cwd = "./dagger/", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("pip install .", cwd = "./dagger/", shell = True, check = True, stdout=stttd, stderr=stttd)
	print("OK!\n")
except:
	print("Failed...\n")


print("Installing fastflood ... ", end = "", flush = True)
try:
	sub.run("rm -r dist", cwd = "./fastflood/", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("rm -r build", cwd = "./fastflood/", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("pip install .", cwd = "./fastflood/", shell = True, check = True, stdout=stttd, stderr=stttd)
	print("OK!\n")
except:
	print("Failed...\n")


print("Installing popscape ... ", end = "", flush = True)
try:
	sub.run("rm -r dist", cwd = "./popscape/", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("rm -r build", cwd = "./popscape/", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("pip install .", cwd = "./popscape/", shell = True, check = True, stdout=stttd, stderr=stttd)
	print("OK!\n")
except:
	print("Failed...\n")
	


print("Done!")
print("\n~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=")
