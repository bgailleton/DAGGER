import matplotlib.pyplot as plt
import numpy as np
import dagger as psc



ny,nx = 256,256
dy,dx = 50,50

Urate = np.zeros((ny,nx)) + 5e-4
Urate[[-1,0],:] = 0
Urate[[0],:] = 1e-4

# topo = np.random.rand(ny,nx)
# topo[-1,:] = 0
# topo[0,:] = 500

dt = 100
Kr = 2e-5
Ks = 5e-5
Kle = 0.1
Kld = 0.1

κ_s = 1.5e-2
κ_r = 1e-2

dep =  4
rshp = (ny,nx)
# Initialising an empty model in the variable ts
ts = psc.trackscape()

# Initialising the topography and its dimensions
ts.init_random(nx, ny,dx,dy,"periodic_EW")
# ts.feed_topo(np.load("topo_ss.npy"))

# ts.init_perlin(nx, ny,dx,dy,"periodic_EW", 5,8, 10,42, False)

otopo = np.copy(ts.get_topo().reshape(rshp))

# FUnctions to set parameters as global homogeneous values (if not initialised, there is a default value)
ts.set_single_Kr(Kr)
ts.set_single_Ks(Ks)
ts.set_single_depcoeff(dep)

ts.set_dt(dt)
ts.set_fluvial_mode(psc.TSC_FLUVIAL.DAVY2009)
# ts.set_fluvial_mode(psc.TSC_FLUVIAL.NONE)
# ts.set_fluvial_mode(psc.TSC_FLUVIAL.FASTSCAPE)
# ts.set_secondary_fluvial_mode(psc.TSC_FLUVIAL.LATERALDAVY)

# ts.set_hillslopes_mode(psc.TSC_HILLSLOPE.CIDRE_NOCRIT)
ts.set_hillslopes_mode(psc.TSC_HILLSLOPE.CIDRE)
ts.set_hillslopes_mode(psc.TSC_HILLSLOPE.HYLANDS)


ts.set_flowtopo_mode(psc.TSC_FLOW_TOPOLOGY.MFD)
# ts.set_flowtopo_mode(psc.TSC_FLOW_TOPOLOGY.SFD)

ts.graph.set_LMR_method(psc.LMR.cordonnier_fill)

ts.set_single_depcoeff(dep)
ts.set_single_Kr(Kr)
ts.set_single_Ks(Ks)

ts.set_single_Kle(Kle)
ts.set_single_Kld(Kld)

ts.set_single_kappa_r(κ_r)
ts.set_single_kappa_s(κ_s)


# Deactivate hilllsopes processes
# ts.hillslopes_off()
ts.strip_sediment()
ts.feed_topo(np.load("topo_ss_hs.npy"))

stopo = ts.get_topo()
ts.set_single_internal_friction(0.5)
ts.set_single_tls(dt*100)

ts.Standalone_hyland_landslides()
node = nx * 200 + 120