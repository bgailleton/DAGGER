import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np
# dem = load_raster("example.tif")
dem = load_raster("../../../../FastFlood2.0/FastFlood2_Boris/data/DEM_1m_crop.tif")

con = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
gf = dag.graph(dem["nx"] * dem["ny"], 8)
gf.init_graph(con)
PPdem = gf.compute_graph("cordonnier_fill", dem['array'].ravel(), con, False, True)
HS = dag.hillshade(con,PPdem)
rshp = [dem['ny'], dem['nx']]
rshp
print(1)

print(2)

SDA = gf.accumulate_constant_downstream_SFD(con, dem['dx'] * dem['dy'])
print(3)


print(4)
gradient = gf.get_links_gradient(con, PPdem)
print(5)

weights = gf.get_link_weights(gradient,0.5)
print(6)
print(weights)

MFA_0 = gf.accumulate_constant_downstream_MFD(con,weights,1)
print(7)

####