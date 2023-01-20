import dagger as dag
import fastflood_dagger as ff
import numpy as np


# number of nodes in the X/Y directions
nx,ny = 200,200
# Resolution in the X/Y directions
dx,dy = 200,200
# generating hte white noise
topo = np.random.rand(ny*nx)+5

con = dag.D8N(nx, ny, dx, dy, 0, 0)
gf = dag.graph(nx * ny, 8)
gf.init_graph(con)


fasft = ff.FF(gf,con,topo.ravel())


