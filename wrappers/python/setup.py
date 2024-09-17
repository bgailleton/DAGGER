import sys, os

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import setuptools
from sys import platform

__version__ = '0.0.1'


class get_pybind_include(object):
	"""Helper class to determine the pybind11 include path

	The purpose of this class is to postpone importing pybind11
	until it is actually installed, so that the ``get_include()``
	method can be invoked. """

	def __init__(self, user=False):
		self.user = user

	def __str__(self):
		import pybind11
		return pybind11.get_include(self.user)


class get_numpy_include(object):
	"""Helper class to determine the numpy include path

	The purpose of this class is to postpone importing numpy
	until it is actually installed, so that the ``get_include()``
	method can be invoked. """

	def __init__(self):
		pass

	def __str__(self):
		import numpy as np
		return np.get_include()



isWindows = None
if platform == "linux" or platform == "linux2":
	isWindows = False
elif platform == "darwin":
	isWindows = False
elif platform == "win32":
	isWindows = True

__version__ = "0.0.8"

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
					include_dirs = ["../../DAGGER", "../../fastflood", "../../popscape",get_pybind_include(), get_pybind_include(user=True), get_numpy_include(),os.path.join(sys.prefix, 'include'), os.path.join(sys.prefix, 'Library', 'include')],
					libraries = LIBBR,
					define_macros = MACROS,
					cxx_std=17,
					language = 'c++',
					extra_compile_args=EXTRA_COMPILE,
					extra_link_args=EXTRA_LINK,
				),
]

class BuildExt(build_ext):
	"""A custom build extension for adding compiler-specific options."""
	c_opts = {
		'msvc': ['/EHsc'],
		'unix': [],
	}

	if sys.platform == 'darwin':
		c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

	def build_extensions(self):
		ct = self.compiler.compiler_type
		opts = self.c_opts.get(ct, [])
		if ct == 'unix':
			opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
			opts.append(cpp_flag(self.compiler))
			if has_flag(self.compiler, '-fvisibility=hidden'):
				opts.append('-fvisibility=hidden')
		elif ct == 'msvc':
			opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
		for ext in self.extensions:
			ext.extra_compile_args = opts
		build_ext.build_extensions(self)
def has_flag(compiler, flagname):
	"""Return a boolean indicating whether a flag name is supported on
	the specified compiler.
	"""
	import tempfile
	with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
		f.write('int main (int argc, char **argv) { return 0; }')
		try:
			compiler.compile([f.name], extra_postargs=[flagname])
		except setuptools.distutils.errors.CompileError:
			return False
	return True


def cpp_flag(compiler):
	"""Return the -std=c++14 compiler flag  and errors when the flag is
	no available.
	"""
	if has_flag(compiler, '-std=c++17'):
		return '-std=c++17'
	else:
		raise RuntimeError('C++17 support is required by xtensor!')

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
	cmdclass={"build_ext": BuildExt},
	zip_safe=False,
	python_requires=">=3.6",
	build_requires=['numpy','pybind11'],
	install_requires=['numpy','pybind11'],
)
