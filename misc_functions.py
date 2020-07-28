# This file contains miscellaneous functions used to create interpolated PSF arrays and simulations
import scipy
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tables import *
import astropy.units as u
import astropy.io.fits as fits
from pathlib import Path
import os
import scipy.ndimage
from scipy.interpolate import RegularGridInterpolator

my_home_path = Path(os.getcwd()) # first, define the paths of each of the files
    
lambda0_m = 575e-9
D = 2.3631
mas_per_lamD = lambda0_m * 360.0 * 3600.0 / (2 * np.pi * D) * 1000    # mas per lambda0/D
as_per_lamD = mas_per_lamD/1000

nzodi = 256
nprop = 256
nipac = 128

pxscl_ipac_lamD = 0.2
pxscl_prop_lamD = 0.1

pxscl_ipac_mas = pxscl_ipac_lamD*mas_per_lamD
pxscl_prop_mas = pxscl_prop_lamD*mas_per_lamD

def closest_PSF(xoff,yoff,
                HLCinterpfun,
                n):
    
    r = np.sqrt(xoff**2+yoff**2)
    theta = np.arctan2(yoff,xoff).to(u.deg)
    
    grid = np.meshgrid(range(n),range(n))
    flattened_grid = np.vstack([grid[0].flatten(),grid[1].flatten()])
    
    pts = np.vstack([ flattened_grid, r.value*np.ones(len(grid[0].flatten())) ]).T
    
    interpped = HLCinterpfun(pts).reshape(n,n).T
    interpped = scipy.ndimage.rotate(interpped,-theta.value,reshape=False,order=3)
    
    return interpped 

def tilt_converter(x,y):
    r = np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x).to(u.deg)
    return r, theta

def mask_im(im,
            pxscl_mas,
            n):
    # Now mask the simulation data
    index = im < im.max()/1e32
    xpix,ypix = np.meshgrid(np.arange(-n/2,n/2),np.arange(-n/2,n/2)) # 256by256 grid
    x = (xpix+.5).flatten()*pxscl_mas
    y = (ypix+.5).flatten()*pxscl_mas
    
    # creating the outer mask
    outer_mask_radius_mas = 9*mas_per_lamD
    index[(np.sqrt((x)**2 + (y)**2)>outer_mask_radius_mas).reshape([n,n])]=True
    
    im_masked = np.ma.masked_array(im,index)
    
    return im_masked

def get_avg_masked(im):
    im_flat = im.flatten()
    max_im = np.max(im)
    count = 1
    total = 0
    for i,val in enumerate(im_flat):
        if np.ma.is_masked(val):
            continue
        count = count + 1
        total += val
    av = total/count
    print('The average of the unmasked data is {:.4f} and the maximum is {:.4f}.'.format(av,max_im))
    return av, max_im

def load_interpped_A(h5fname):
    
    interpped_path = '/groups/douglase/Interpped_PSFs/' # directory containing the h5 files
    fpath = interpped_path + h5fname
    print('Loading interpolated PSFs from ' + fpath + '.')

    start = time.time() # start time for loading in the array
    
    h5f = open_file(fpath, mode="r")
    tab = h5f.root.interpolated_library # points to the table containing the data

    shape = (tab[0]['array'].shape[0],tab.shape[0])
    A = np.array(np.zeros(shape)) # this might take a lot of memory!

    #populate an array
    for i in range(shape[1]):
        A[:,i] = np.float_(tab[i]['array'])
    end = time.time()
    h5f.close()
    
    print('Matrix loaded in {:.4f}'.format(end-start))
    return A

def run_sim_forloop(zodi,fname,n):
    zodi_flat = zodi.flatten()
    
    interpped_path = '/groups/douglase/Interpped_PSFs/'
    fpath = interpped_path + fname
    print('Loading interpolated PSFs from ' + fpath + '.')
    
    h5f = open_file(fpath, mode="r")
    tab = h5f.root.interpolated_library # points to the table containing the data
    
    start=time.time()
    sim = np.zeros(n**2)
    for i,zodi_val in enumerate(zodi_flat):
        sim += tab[i]['array']*zodi_val
    sim = sim.reshape(n,n)
    end = time.time()
    h5f.close()
    
    print('Simulation done in {:.4f}'.format(end-start))
    return sim
        
    
def run_sim_standard(zodi,HLCinterpfun,n):
    zodi_flat = zodi.flatten()
    
    zodi_pixscale = 3.529*u.milliarcsecond
    xpix,ypix = np.meshgrid(np.arange(-nzodi/2,nzodi/2),np.arange(-nzodi/2,nzodi/2)) # 256by256 grid
    x = (xpix+.5).flatten()*zodi_pixscale
    y = (ypix+.5).flatten()*zodi_pixscale
    
    print('Running simulation...')
    start = time.time()
    sim = np.zeros(shape=(n,n))
    for i,zodi_val in enumerate(zodi_flat):
        interpped = closest_PSF(x[i],y[i],HLCinterpfun,n)
        sim += zodi_val*interpped
    end = time.time()
    
    print('Simulation done in {:.4f}'.format(end-start))
    return sim


def setup_ipac_interpfun(method='nearest'): # this functions automaticall sets up all the interpolating functions
    # start with the ipac psfs
    ipac_psfs_file_path = my_home_path/'offset_psfs'/'IPAC_HLC_PSFs'/'20180718_hlc_nfov_PSFs_1Dcube.fits'
    ipac_psfs = fits.getdata(ipac_psfs_file_path)
    ipac_psfs = np.moveaxis(ipac_psfs,0,-1) # shift the axis of the datacube such that the final entry is the psf number

    # get the offset sampling of the ipac psfs from the info file
    ipac_info_file_path = my_home_path/'offset_psfs'/'IPAC_HLC_PSFs'/'20180718_hlc_nfov_PSFs_1Dcube_info.fits'
    ipac_offsets = fits.getdata(ipac_info_file_path)
    offsets_lamD_ipac = ipac_offsets[:,1]
    offsets_mas_ipac = offsets_lamD_ipac*mas_per_lamD

    # setup the interpolating functions
    x = range(nipac)
    y = range(nipac)
    z = offsets_mas_ipac
    ipac_interpfun = RegularGridInterpolator((x, y, z), 
                                             ipac_psfs,
                                             bounds_error = False,
                                             method = method,
                                             fill_value = 0)
    return ipac_interpfun
    
    
def setup_prop_interpfun(method='nearest'):
    # next use the  the non DM proper psfs
    prop_psfs_file_path = my_home_path/'offset_psfs'/'PROPER_HLC_PSFs'/'20200721_HLC_1D_offaxis_psfs.fits'
    prop_psfs = fits.getdata(prop_psfs_file_path)
    prop_psfs = np.moveaxis(prop_psfs,0,-1)

    # make the offset sampling of the PROPER PSFs
    start = 0
    stop = 11
    step = 0.05
    offsets_lamD_prop = np.arange(start,stop+step,step)
    offsets_mas_prop = offsets_lamD_prop*mas_per_lamD

    # setup the interpolating functions
    x = range(nprop)
    y = range(nprop)
    z = offsets_mas_prop
    prop_interpfun = RegularGridInterpolator((x, y, z), 
                                               prop_psfs,
                                               bounds_error = False,
                                               method = method,
                                               fill_value = 0)

    return  prop_interpfun
    
def setup_prop_dm_interpfun(method='nearest'):
    # next use the proper psfs with DMs 
    prop_DM_psfs_file_path = my_home_path/'offset_psfs'/'PROPER_HLC_PSFs'/'20200720_HLC_1D_offaxis_psfs_DMs.fits'
    prop_DM_psfs = fits.getdata(prop_DM_psfs_file_path)
    prop_DM_psfs = np.moveaxis(prop_DM_psfs,0,-1)

    # make the offset sampling of the PROPER PSFs
    start = 0
    stop = 11
    step = 0.05
    offsets_lamD_prop = np.arange(start,stop+step,step)
    offsets_mas_prop = offsets_lamD_prop*mas_per_lamD

    # setup the interpolating functions
    x = range(nprop)
    y = range(nprop)
    z = offsets_mas_prop
    prop_dm_interpfun = RegularGridInterpolator((x, y, z), 
                                                prop_DM_psfs,
                                                bounds_error = False,
                                                method = method,
                                                fill_value = 0)

    return prop_dm_interpfun


def apply_det(image,
              ps_det=0.0211*u.arcsec/u.pixel,#CGI_Imaging_Pixel_Scale
              ps_input=0.005*u.arcsec/u.pixel,#OS6 PSF simulations
              det_shape=(150,150),
              verbose=False):
    if hasattr(image, 'unit'):
        unit=image.unit
    else:
        unit=1
    #zoomed= resize(image/image.sum(),shape,anti_aliasing=True,order=1) #(image, output_shape, order=1, mode='reflect', cval=0, clip=True, preserve_range=False, anti_aliasing=True, anti_aliasing_sigma=None)
    #return 
    zoom = (ps_input/ps_det).decompose().value
    resampled_image = scipy.ndimage.interpolation.zoom(image, 
                                                       zoom,
                                                       output=image.dtype,
                                                       order=3)
    
    lx, ly = resampled_image.shape
    if verbose:
        print("zoom: "+str(zoom))
    # crop down to match size of detector:
    if det_shape == None:
        lx_w,ly_w = resampled_image.shape
    else:
        lx_w, ly_w = det_shape
    border_x = np.abs(lx - lx_w) // 2
    border_y = np.abs(ly - ly_w) // 2


    if zoom<1:
        new_im = np.zeros([lx_w, ly_w])
        new_im[border_x:border_x + resampled_image.shape[0],
                                        border_y:border_y + resampled_image.shape[1]] = resampled_image
    else:
        new_im=resampled_image[border_x:border_x + lx_w, border_y:border_y + ly_w]
    
    renormed_im = new_im*image.sum()/new_im.sum()
    #print(renormed_im.sum())
    return resampled_image, renormed_im #check this conserves flux well enough
