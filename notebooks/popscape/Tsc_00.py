import dagger as psc
import matplotlib.pyplot as plt
import numpy as np


ny,nx = 256,256
dy,dx = 200,200

Urate = np.zeros((ny,nx)) + 5e-4
Urate[round(2*ny/3):,:] = 0

dt = 1000
Kr = 1e-5
Ks = 2e-5
dep = 1000
rshp = (ny,nx)
# Initialising an empty model in the variable ts
ts = psc.trackscape()
# Initialising the topography and its dimensions
ts.init_random(nx, ny,dx,dy,"periodic_EW")
# ts.init_perlin(nx, ny,dx,dy,"periodic_EW", 5,8, 200,42, False)

# FUnctions to set parameters as global homogeneous values (if not initialised, there is a default value)
ts.set_single_Kr(Kr)
ts.set_single_Ks(Ks)
ts.set_single_depcoeff(dep)

ts.set_dt(dt)
ts.set_fluvial_mode(psc.TSC_FLUVIAL.DAVY2009)
# ts.set_fluvial_mode(psc.TSC_FLUVIAL.FASTSCAPE)
ts.set_flowtopo_mode(psc.TSC_FLOW_TOPOLOGY.SFD)

ts.graph.set_LMR_method(psc.LMR.cordonnier_carve)
ndt = 100000
nupdate = 100

# Main loop
for i in range(ndt):
    if(i<500):
        ts.set_dt(dt/10)
    else:
        ts.set_dt(dt)     

    ts.run()
    # ts.run_SFD(dt)
    ts.external_uplift(Urate,dt, False)