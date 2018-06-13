import numpy as np
import xarray as xr

def rotate_3D(ds, u_unrot_name='', v_unrot_name='', w_unrot_name=''):
    '''
    Function performing a 3D rotation of u,v,w from sonic anemometer
    performing the first two rotation, i.e. (1) aligning along and cross wind
    components and (2) mean vertical wind equals zero

    Developed and written by Christoph Thomas
    Dept. of Forest Science, Oregon State University, 2006
    Converted to python by Karl lapo
    Micrometeorology Group, University of Bayreuth, 2018
    '''

    # Calculate averages
    u_mean_unrot = np.nanmean(ds[u_unrot_name])
    v_mean_unrot = np.nanmean(ds[v_unrot_name])
    w_mean_unrot = np.nanmean(ds[w_unrot_name])

    # Unpack unrotated wind vectors
    u_unrot = ds[u_unrot_name]
    v_unrot = ds[v_unrot_name]
    w_unrot = ds[w_unrot_name]

    # Performing the first rotation (around the w-axis)
    alpha_rot = np.arctan(v_mean_unrot / u_mean_unrot)
    if u_mean_unrot < 0:
        alpha_rot = alpha_rot + np.pi
    elif v_mean_unrot < 0:
        alpha_rot = alpha_rot + 2 * np.pi

    u_rot = (u_unrot * np.cos(alpha_rot)) + (v_unrot * np.sin(alpha_rot))
    v_rot = (-u_unrot * np.sin(alpha_rot)) + (v_unrot * np.cos(alpha_rot))

    # Calculating the second rotation angle beta (mean w equals zero)
    beta_rot = np.arctan2(w_mean_unrot, np.nanmean(u_rot))

    # Performing the second rotation (around the v-axis)
    u_rot = (u_rot * np.cos(beta_rot))  + (w_unrot * np.sin(beta_rot))
    w_rot = -(u_rot * np.sin(beta_rot)) + (w_unrot * np.cos(beta_rot))

    # List of angles
    alpha_rot = np.rad2deg(alpha_rot)
    beta_rot = np.rad2deg(beta_rot)

    # Pack back into the xarray dataset
    ds[u_unrot_name + '_rotated'] = u_rot
    ds[v_unrot_name + '_rotated'] = v_rot
    ds[w_unrot_name + '_rotated'] = w_rot
    ds['alpha'] = beta_rot
    ds['beta'] = alpha_rot

    return(ds)
