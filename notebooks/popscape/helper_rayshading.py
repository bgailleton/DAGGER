import tempfile

import matplotlib.pyplot as plt
import numpy as np
import rasterio as rio
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import rpy2.robjects.packages as rpackages

rpy2.robjects.numpy2ri.activate()

import scabbard as scb
from scipy.ndimage import gaussian_filter


def rayshade(tackmod, ny, nx, ngauss_topo = 3, ngauss_water = 1, threshold_A = 1e4, img_path=None, zscale=10, fov=0, theta=135, zoom=0.75, phi=45, windowsize=(1000, 1000), shader = "desert"):
    
    # Output path.
    if not img_path:
        img_path = tempfile.NamedTemporaryFile(suffix='.png').name
    
    # Import needed packages.
    rayshader = rpackages.importr('rayshader')
    
    z = tackmod.get_topo().reshape(ny,nx)
    z = gaussian_filter(z,ngauss_topo)
    Qw = tackmod.get_Qw().reshape(ny,nx)
    if(ngauss_water > 0):
        Qw = gaussian_filter(Qw,ngauss_water)
    
    water = (Qw > threshold_A)
    
    
    # Convert array to matrix.
    z = np.asarray(z)
    rows, cols = z.shape
    z_mat = ro.r.matrix(z, nrow=rows, ncol=cols)
    ro.globalenv['elmat'] = z_mat
    
    water = np.asarray(water)[::-1]
    water_mat = ro.r.matrix(water, nrow=rows, ncol=cols)
    ro.globalenv['wamat'] = water_mat
    
    # Save python state to r.
    ro.globalenv['img_path'] = img_path
    ro.globalenv['zscale'] = zscale
    ro.globalenv['fov'] = fov
    ro.globalenv['theta'] = theta
    ro.globalenv['zoom'] = zoom
    ro.globalenv['phi'] = phi
    ro.globalenv['windowsize'] = ro.IntVector(windowsize)
    
    # Do the render.
    ro.r(f'''
        elmat %>%
          sphere_shade(texture = "{shader}") %>%
          #sphere_shade(texture = "imhof1") %>%
          #sphere_shade(texture = "imhof2") %>%
          #sphere_shade(texture = "imhof3") %>%
          #sphere_shade(texture = "imhof4") %>%
          #sphere_shade(texture = "bw") %>%
          #sphere_shade(texture = "unicorn") %>%
          add_water(wamat, color = "{shader}") %>%
          add_shadow(ray_shade(elmat, zscale = 1), 0.7) %>%
          add_shadow(ambient_shade(elmat), 0) %>%
          plot_3d(elmat, zscale = zscale, fov = fov, theta = theta, zoom = zoom, phi = phi, windowsize = windowsize)
        Sys.sleep(0.2)
        render_snapshot(img_path)
    ''')
    
    # Return path.
    return img_path