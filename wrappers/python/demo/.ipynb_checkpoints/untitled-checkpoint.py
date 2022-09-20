import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np

dem = load_raster("example.tif")

con = dag.D4N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
gf = dag.graphD4(dem["nx"] * dem["ny"], 4)
gf.init_graph(con)
gf.set_LMR_method(dag.LMR.cordonnier_fill)
for i in range(50):
    print(i)
    PPdem = gf.compute_graph(dem['array'].ravel(), con, False, True)
# HS = dag.hillshade(con,PPdem)
rshp = [dem['ny'], dem['nx']]