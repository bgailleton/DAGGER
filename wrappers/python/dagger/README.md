# `python` bindings to `DAGGER`

In this folder, you can find the `python` bindings to `DAGGER` using `pybind11`. A `conda-forge` package will be available in the future, but so far you can easily install it in a `conda` environment.

## Quickstart if comfortable with `conda`

1) clone the full repository

2) `conda install -c conda-forge rasterio gdal scipy numpy matplotlib jupyterlab numpy pybind11 cmake scikit-build ipympl numba`

3) `python setup.py install` (if you do not have a compiler installed, `conda install -c conda-forge clang` should solve the issue after a terminal reboot)

4) cd to `demo` folder and run `jupyter lab` to open a jupyter workspace and explore the example notebooks


## I am beginning with `python/conda/...`

I'll write more detailed instructions later but long story short `python` requires an **environment** to know which packages are available. To simplify its management, `conda` creates "boxes" inside the computer isolating a `python` environment in order to install anything you need in it without worrying about compatibility issues. It also provide easy ways to install and distribute libraries (not only in python but also R/Julia/C++/Whatever).

Multiple implementations exists, but to guarantee full FOSS I recommend to install [`mambaforge`](https://github.com/conda-forge/miniforge#download). Once install, you can run in your terminal the `conda` commands.

First let's create a new environment (i.e. creating the box, only needed once):

`conda create -n daggerenv`

Then we need to activate it (i.e. enter the box, needed every time you wanna use the environment inside it):

`conda activate daggerenv`

Now you can go to the section above to complete the install.
