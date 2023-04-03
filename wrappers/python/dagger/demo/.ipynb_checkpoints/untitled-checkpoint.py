import dagger as dag
from helper import load_raster
import matplotlib.pyplot as plt
import numpy as np
import scabbard as scb


gf.set_LMR_method(dag.LMR.dagger_carve )
PPdem = gf.compute_graph(dem['array'].ravel(), False, False)
HS = dag.hillshade(con,PPdem)
rshp