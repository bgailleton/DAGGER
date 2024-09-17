# Demo notebooks and script using DAGGER in `python`

In this folder, you will find demo scripts and jupyter notebook using different aspects of `DAGGER` in `python`. This is the language I spent the most time binding as one of the most mature and used in Geosciences (And eventually because this is my everyday language with `c++`).

In order to run these examples out of the box, I reccomend you install the folowing packages, they are not mandatory to use `DAGGER`, but will help make these example nice and easy (e.g. loading tif DEMs, running stats and general anaylis, visualisation, data I/O ...).


```
# If you have conda:
conda install -c conda-forge matplotlib ipympl jupyter-lab rasterio gdal numba

# if you directly installed mambaforge:
mamba install matplotlib ipympl jupyter-lab rasterio gdal numba
```

## Notebooks and scripts

- `drainage_area`: shows how to calculate drainage area using `DAGGER`. It illustrates the differences between multiple flow direction (MFD) and single flow direction (SFD); how to preprocess a DEM and why it is important;
