import sys, os
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension
from sys import platform

__version__ = '0.0.8'

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path."""
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

class get_numpy_include(object):
    """Helper class to determine the numpy include path."""
    def __str__(self):
        import numpy as np
        return np.get_include()

isWindows = platform == "win32"

# Compilation flags and macros
MACROS = [('VERSION_INFO', __version__), ("DAGGER_FT_PYTHON", None), ("BOOST_AVAILABLE", None)]
EXTRA_COMPILE = ['/Ox'] if isWindows else ['-O3', '-Wall', '-std=c++17']
EXTRA_LINK = [] if isWindows else ['-O3']
LIBBR = []

if "--exp" in sys.argv:
    MACROS.append(('OPENMP_YOLO', None))
    EXTRA_COMPILE.append('-fopenmp')
    sys.argv.remove("--exp")
    LIBBR.append("gomp")

ext_modules = [
    Pybind11Extension(
        "dagger",
        ["main.cpp"],
        include_dirs=[
            "../../DAGGER", "../../fastflood", "../../popscape",
            get_pybind_include(), get_pybind_include(user=True), get_numpy_include()
        ],
        libraries=LIBBR,
        define_macros=MACROS,
        cxx_std=17,
        language='c++',
        extra_compile_args=EXTRA_COMPILE,
        extra_link_args=EXTRA_LINK,
    ),
]

class BuildExt(build_ext):
    """Custom build extension to add compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '/std:c++17'],  # Ensures C++17 is used in MSVC
        'unix': ['-std=c++17'],
    }

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

def has_flag(compiler, flagname):
    """Check if the compiler supports a specific flag."""
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main() { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True

setup(
    name="DAGGER",
    version=__version__,
    author="Boris Gailleton",
    author_email="boris.gailleton@univ-rennes.fr",
    url="https://github.com/bgailleton/DAGGER",
    description="Directed Acyclic Graph for digital topoGrahic manipulations Eventually cRoss-platform",
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=['numpy', 'pybind11'],
)
