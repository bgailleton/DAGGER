import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np
import scabbard as scb

nx = 5
ny = 10
dx,dy = 1,1

S0 = 1

BC = np.ones((ny,nx), dtype=np.uint8)
BC[0,:] = 0
BC[:,[0,-1]] = 0
BC[-1,1:9] = 4

grid = scb.slope_RGrid(nx,ny,dx,dy, slope = S0, noise_magnitude=0., EW = "out", S = "out", N = "noflow")

grid.Z2D[5,1] = 0
grid.Z2D[5,3] = 0
grid.Z2D[2:7,1:4] -= 5

grid.compute_graphcon()
grid.con.set_custom_boundaries(BC.ravel())

grid.graph.set_LMR_method(dag.LMR.dagger_carve )
PPdem = grid.graph.compute_graph(grid.Z, False, False)

A = grid.graph.accumulate_constant_downstream_SFD(1)

fig,ax = plt.subplots()
ax.imshow(A.reshape(grid.rshp))
fig.show()
B = np.copy(A)
B[A>=1] =1
print(np.sum(B))

input()

