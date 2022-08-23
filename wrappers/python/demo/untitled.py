import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np
import time
# dem = load_raster("example.tif")
print(1)
dem = load_raster("../../../../FastFlood2.0/FastFlood2_Boris/data/DEM_1m_crop.tif")
print(2)

con = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
print(3)
gf = dag.graph(dem["nx"] * dem["ny"], 8)
print(4)
gf.init_graph(con)
# for i in range(50):

print(5)
PPdem = gf.compute_graph("cordonnier_fill", dem['array'].ravel(), con, False, True)
# time.sleep(1)
print(6)
HS = dag.hillshade(con,PPdem)
print(7)
rshp = [dem['ny'], dem['nx']]
rshp
print(11)

print(12)

SDA = gf.accumulate_constant_downstream_SFD(con, dem['dx'] * dem['dy'])
print(13)


print(14)
gradient = gf.get_links_gradient(con, PPdem)
print(15)

weights = gf.get_link_weights(gradient,0)
print(16)
# print(weights)
st = time.time()

# for i in range(50):
MFA_0 = gf.accumulate_constant_downstream_MFD(con,weights,900)
np.save("temp.npy", MFA_0)
print("took::", (time.time() - st)/50)

print(np.unique(MFA_0))
print(17)

####