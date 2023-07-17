.. _quickstart:

Quick start
===========

``c++``
--------

As a header-only ``c++`` library, you simply need to include the right headers in your own ``c++`` code to use ``DAGGER`` and include the files in your local compilation system. The name of the header files depends on the modules used and can be found in the usage part of the compilation.

``python``
----------

By far the easiest way to install the ``python`` package (direct bindings) is by using a ``conda`` environment and ``conda-forge``. The package is (erm... will be) installable with ``conda install -c conda-forge dagger``. Guidance on Anaconda or package building from source can be found in the `installation`_ section.

You can then find numerous ``jupyter notebooks`` illustrating the use of ``DAGGER``  in the `demo folder <https://github.com/bgailleton/DAGGER/tree/main/wrappers/python/demo>`_.

``Julia``
-----------

There is no pre-compiled package for using ``DAGGER`` with Julia (yet). Instructions on how to use it can be found `there <https://github.com/bgailleton/DAGGER/tree/main/wrappers/julia>`_, and ``demo`` `here <https://github.com/bgailleton/DAGGER/tree/main/wrappers/julia/demo>`_.

``R``
----------

TODO

``JS/webassembly``
-------------------

TODO

``MATLAB``
-----------

MATLAB's structure and nature makes the binding with ``c++`` difficult and inefficient - at least without the expensive compiler toolbox. A proof of concept showing it is possible can be found there (TODO) and MATLAB bindings will be limited to specific projects.
