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


# none = True if (("-g" in sys.argv) == False and ("-d" in sys.argv) == False and ("-p" in sys.argv) == False ) else True
# dagger = True if ("-d" in sys.argv or none) else False
# graphflood = True if ("-g" in sys.argv or none) else False
# popscape = True if ("-p" in sys.argv or none) else False

print(f"\n\n{bcolors.YELLOW}~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~={bcolors.ENDC}")
print(f"{bcolors.HEADER}I am going to install/update dagger suite codes from source in python.{bcolors.ENDC}\n")
print(f"{bcolors.BOLD}DAGGER core:{bcolors.ENDC} {bcolors.OKGREEN }{'included'}{bcolors.ENDC} ")
print(f"{bcolors.BOLD}graphflood core:{bcolors.ENDC} {bcolors.OKGREEN }{'included'}{bcolors.ENDC} ")
print(f"{bcolors.BOLD}popscape core:{bcolors.ENDC} {bcolors.OKGREEN}{'included'}{bcolors.ENDC} ")


print("\n")

if("-v" in sys.argv):
	stttd = None
else:
	print(f"{bcolors.WARNING}use option -v to get the whole install log{bcolors.ENDC}\n")
	stttd = sub.DEVNULL

print(f"{bcolors.BOLD}Installing ... {bcolors.ENDC}", end = "", flush = True)
run_test = True
try:
	sub.run("pip uninstall -y dagger", shell = True, check = False, stdout = stttd, stderr = stttd)
	sub.run("pip uninstall -y popscape-dagger", shell = True, check = False, stdout = stttd, stderr = stttd)
	sub.run("pip uninstall -y fastflood-dagger", shell = True, check = False, stdout = stttd, stderr = stttd)
	sub.run("rm -r dist", cwd = "./", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("rm -r build", cwd = "./", shell = True, check = False, stdout=stttd, stderr=stttd)
	sub.run("pip install .", cwd = "./", shell = True, check = True, stdout=stttd, stderr=stttd)
	print(f"{bcolors.OKGREEN}OK!{bcolors.ENDC}\n")
except:
	print(f"{bcolors.FAIL}Failed...{bcolors.ENDC}\n")
	run_test = False

if(run_test):
	print(f"{bcolors.BOLD}Running tests ... {bcolors.ENDC}", end = "", flush = True)
	print(f"{bcolors.YELLOW}TODO{bcolors.ENDC}\n")


print(f"{bcolors.BOLD}Done!{bcolors.ENDC}")
print(f"\n{bcolors.YELLOW}~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~={bcolors.ENDC}")
