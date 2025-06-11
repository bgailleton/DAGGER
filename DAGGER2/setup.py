from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11.setup_helpers import ParallelCompile
import pybind11

# Optional: Enable parallel compilation
ParallelCompile("NPY_NUM_BUILD_JOBS").install()

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "dg2",  # Module name
        sources=[
            "main.cpp",  # Your main cpp file
        ],
        include_dirs=[
            # Path to pybind11 headers
            pybind11.get_cmake_dir() + "/../../../include",
            ".",  # Current directory for dg2_array.hpp
        ],
        language='c++',
        cxx_std=17,  # C++14 standard (required for pybind11),
        define_macros = [("DG2_WITH_PYTHON", None)]
    ),
]

setup(
    name="dg2",
    version="0.0.1",
    author="Boris Gailleton",
    author_email="boris.gailleton@univ-rennes.fr",
    description="Efficient NumPy-C++ array interoperability utilities",
    long_description="""
    WIP
    """,
    long_description_content_type="text/plain",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.15",
    ],
    extras_require={
        "test": ["pytest", "numpy"],
        "dev": ["pytest", "numpy", "build", "twine"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords="numpy, c++, arrays, performance, interop",
    project_urls={
        "Source": "https://github.com/yourusername/dg2-arrays",
        "Bug Reports": "https://github.com/yourusername/dg2-arrays/issues",
    },
)
