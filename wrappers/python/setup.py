import sys, os

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


from sys import platform

isWindows = None
if platform == "linux" or platform == "linux2":
	isWindows = False
elif platform == "darwin":
	isWindows = False
elif platform == "win32":
	isWindows = True

__version__ = "0.0."

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#	 say from a submodule.
#
# Note:
#	 Sort input source files if you glob sources to ensure bit-for-bit
#	 reproducible builds (https://github.com/pybind/python_example/pull/53)
MACROS = [('VERSION_INFO', __version__), ("DAGGER_FT_PYTHON", None), ("BOOST_AVAILABLE", None)]

if(isWindows):
	EXTRA_COMPILE = ['/Ox']
	EXTRA_LINK = []

else:
	EXTRA_COMPILE = ['-O3', '-Wall']
	EXTRA_LINK = ['-O3']


# LIBBR = ["boost_system"]
LIBBR = []


if "--exp" in sys.argv:
	MACROS.append(('OPENMP_YOLO', None))
	EXTRA_COMPILE.append('-fopenmp')
	EXTRA_COMPILE.append('-foffload=nvptx-none')
	sys.argv.remove("--exp")
	import os
	os.environ["CC"] = "g++-12"
	LIBBR.append("gomp")

ext_modules = [
		Pybind11Extension(
					"dagger",
					["main.cpp"],
					include_dirs = ["../../DAGGER", "../../fastflood", "../../popscape"],
					libraries = LIBBR,
					define_macros = MACROS,
					cxx_std=17,

					 extra_compile_args=EXTRA_COMPILE,
					 extra_link_args=EXTRA_LINK,

				),
]

setup(
		name="DAGGER",
		version=__version__,
		author="Boris Gailleton",
		author_email="boris.gailleton@univ-rennes.fr",
		url="https://github.com/bgailleton/DAGGER",
		description="Directed Acyclic Graph for digital topoGrahic manipulations Eventually cRoss-platform",
		long_description="",
		ext_modules=ext_modules,
		extras_require={"test": "pytest"},
		# Currently, build_ext only provides an optional "highest supported C++
		# level" feature, but in the future it may provide more features.
		cmdclass={"build_ext": build_ext},
		zip_safe=False,
		python_requires=">=3.6",
)
