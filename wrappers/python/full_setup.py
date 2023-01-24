'''
This script runs a full install of the dagger-related codes
B.G.
'''

import subprocess as sub
import sys

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
	YELLOW = '\033[33;49m'


none = True if (("-g" in sys.argv) == False and ("-d" in sys.argv) == False and ("-p" in sys.argv) == False ) else True
dagger = True if ("-d" in sys.argv or none) else False
graphflood = True if ("-g" in sys.argv or none) else False
popscape = True if ("-p" in sys.argv or none) else False

print(f"\n\n{bcolors.YELLOW}~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~={bcolors.ENDC}")
print(f"{bcolors.HEADER}I am going to install/update dagger suite codes from source in python.{bcolors.ENDC}\n")
print(f"{bcolors.BOLD}DAGGER core:{bcolors.ENDC} {bcolors.OKGREEN if dagger else bcolors.FAIL}{'yes' if dagger else 'no'}{bcolors.ENDC} ")
print(f"{bcolors.BOLD}graphflood core:{bcolors.ENDC} {bcolors.OKGREEN if dagger else bcolors.FAIL}{'yes' if graphflood else 'no'}{bcolors.ENDC} ")
print(f"{bcolors.BOLD}popscape core:{bcolors.ENDC} {bcolors.OKGREEN if popscape else bcolors.FAIL}{'yes' if dagger else 'no'}{bcolors.ENDC} ")


print("\n")

if("-v" in sys.argv):
	stttd = None
else:
	print(f"{bcolors.WARNING}use option -v to get the whole install log{bcolors.ENDC}\n")
	stttd = sub.DEVNULL

if(dagger or none):
	print(f"{bcolors.BOLD}Installing dagger ... {bcolors.ENDC}", end = "", flush = True)
	try:
		sub.run("pip uninstall -y dagger", shell = True, check = False, stdout = stttd, stderr = stttd)
		sub.run("rm -r dist", cwd = "./dagger/", shell = True, check = False, stdout=stttd, stderr=stttd)
		sub.run("rm -r build", cwd = "./dagger/", shell = True, check = False, stdout=stttd, stderr=stttd)
		sub.run("pip install .", cwd = "./dagger/", shell = True, check = True, stdout=stttd, stderr=stttd)
		print(f"{bcolors.OKGREEN}OK!{bcolors.ENDC}\n")
	except:
		print(f"{bcolors.FAIL}Failed...{bcolors.ENDC}\n")

if(graphflood or none):
	print(f"{bcolors.BOLD}Installing graphflood ... {bcolors.ENDC}", end = "", flush = True)
	try:
		sub.run("pip uninstall -y fastflood-dagger", shell = True, check = False, stdout = stttd, stderr = stttd)
		sub.run("rm -r dist", cwd = "./fastflood/", shell = True, check = False, stdout=stttd, stderr=stttd)
		sub.run("rm -r build", cwd = "./fastflood/", shell = True, check = False, stdout=stttd, stderr=stttd)
		sub.run("pip install .", cwd = "./fastflood/", shell = True, check = True, stdout=stttd, stderr=stttd)
		print(f"{bcolors.OKGREEN}OK!{bcolors.ENDC}\n")
	except:
		print(f"{bcolors.FAIL}Failed...{bcolors.ENDC}\n")

if(popscape or none):
	print(f"{bcolors.BOLD}Installing popscape ... {bcolors.ENDC}", end = "", flush = True)
	try:
		sub.run("pip uninstall -y popscape-dagger", shell = True, check = False, stdout = stttd, stderr = stttd)
		sub.run("rm -r dist", cwd = "./popscape/", shell = True, check = False, stdout=stttd, stderr=stttd)
		sub.run("rm -r build", cwd = "./popscape/", shell = True, check = False, stdout=stttd, stderr=stttd)
		sub.run("pip install .", cwd = "./popscape/", shell = True, check = True, stdout=stttd, stderr=stttd)
		print(f"{bcolors.OKGREEN}OK!{bcolors.ENDC}\n")
	except:
		print(f"{bcolors.FAIL}Failed...{bcolors.ENDC}\n")
	


print(f"{bcolors.BOLD}Done!{bcolors.ENDC}")
print(f"\n{bcolors.YELLOW}~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~={bcolors.ENDC}")
