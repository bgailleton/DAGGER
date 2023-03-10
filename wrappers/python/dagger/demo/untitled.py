import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np

dem = load_raster("putna_FAKE_20.tif")
dem['array'] = np.random.rand(dem['nx'] * dem['ny'])

con = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
gf = dag.graph(con)
gf.init_graph()

gf.set_LMR_method(dag.LMR.cordonnier_carve)
PPdem = gf.compute_graph(dem['array'].ravel(), True, True)