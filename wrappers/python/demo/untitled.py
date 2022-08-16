import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np

dem = load_raster("putna_50_NDF.tif")
dem['array'] += np.random.rand(dem["ny"], dem["nx"]) * 100

con = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
gf = dag.graph(dem["nx"] * dem["ny"], 8)

gf.init_graph(con)


PPdem = gf.compute_graph("cordonnier_fill", dem['array'].ravel(), con, False, True)

gradient = gf.get_links_gradient(con, PPdem)
weights = gf.get_link_weights(gradient,1)
A = gf.accumulate_constant_downstream_MFD(con,weights,900)

print(np.unique(weights))







####