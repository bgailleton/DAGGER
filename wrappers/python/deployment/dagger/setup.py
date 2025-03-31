#!/usr/bin/env python

"""The setup script."""

# pybind11 is used for c++ integration
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages
import os
import sys
import platform

__version__ = "0.0.17"



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





#############################################
# Admin inputs
#############################################

# including readme info
with open('README.rst') as readme_file:
    readme = readme_file.read()

# including historical changes and additions
with open('HISTORY.rst') as history_file:
    history = history_file.read()

#############################################
# Requirements
#############################################

# Need click for common CLI (although there are not any yet)
# Need pybind11 for compilation
requirements = ['Click>=7.0', 'pybind11', 'numpy>=2']

# Testing happens internally before deployment
test_requirements = [ ]

#############################################
# C++ compiler flag, libraries, speciations and all
#############################################

# Get the operating system name as a string
os_name = platform.system()

# Macros to enable python wrappers for I/O (may become platform dependent for opti or experimental features?)
MACROS = [('VERSION_INFO', __version__), ("DAGGER_FT_PYTHON", None)]

# Platform specific flags
if(os_name == 'Windows'):
    EXTRA_COMPILE = ['/Ox', '/EHsc', '/std:c++17']
    EXTRA_LINK = ['/Ox']
else:
    EXTRA_COMPILE = ['-O3', '-Wall', '-Wextra', '-std=c++17']
    EXTRA_LINK = ['-O3']

# Extra libaries to fetch. I'd like to keep that minimal but we never know
LIBBR = []

# Actual compilation
ext_modules = [
        Pybind11Extension(
                    "dagger",
                    ["main.cpp"], # DAGGER is header-only, so only the main cpp file containing definitions has to be compiled
                    include_dirs = ["includes", get_pybind_include(), get_pybind_include(user=True), get_numpy_include(),os.path.join(sys.prefix, 'include'), os.path.join(sys.prefix, 'Library', 'include')], # Deployment process fetch all the required headers and gather them in `includes`
                    libraries = LIBBR, # Additional libraries
                    define_macros = MACROS, # -D stuff
                    cxx_std=17, # I need c++17 standard for most of my tools
                    extra_compile_args=EXTRA_COMPILE, # compiler flags
                    extra_link_args=EXTRA_LINK, # compiler flags
                ),
]

#############################################
# Python setup
#############################################

## Main setup function
setup(
    author="Boris Gailleton",
    author_email='boris.gailleton@univ-rennes.fr',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    description="DAG tools to process numerical topography and landscape evolution models",
    entry_points={
        'console_scripts': [
            'dagger=dagger.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='dagger',
    name='daggerpy',
    test_suite='tests',
    tests_require=test_requirements,
    ext_modules=ext_modules,
    url='https://github.com/bgailleton/DAGGER',
    version=__version__,
    zip_safe=False,
    cmdclass={"build_ext": build_ext},
)
