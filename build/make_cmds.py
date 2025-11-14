#! /usr/bin/env python

"""Makes CMDs of GST measurements.

Authors
-------
    Meredith Durbin, February 2018

Use
---
    This script is intended to be executed from the command line as
    such:
    ::
        python make_cmds.py ['filebase']
    
    Parameters:
    (Required) [filebase] - Path to .phot.hdf5 file
"""
import matplotlib
matplotlib.use('Agg') # check this
import matplotlib.pyplot as plt
import numpy as np
import wpipe as wp
import dask.dataframe as dd
import os
import pandas as pd
import traceback
from astropy.io import fits
import argparse
from astropy.wcs import WCS
from pathlib import Path
import time


def register(task):
    _temp = task.mask(source="*", name="start", value=task.name)
    _temp = task.mask(source="*", name="hdf5_ready", value="*")

# global photometry values
# first 11 columns in raw dolphot output
global_columns = ['ext','chip','x','y','chi_gl','snr_gl','sharp_gl', 
                  'round_gl','majax_gl','crowd_gl','objtype_gl']

# dictionary mapping text in .columns file to column suffix
colname_mappings = {
    'counts,'                            : 'count',
    'sky level,'                         : 'sky',
    'Normalized count rate,'             : 'rate',
    'Normalized count rate uncertainty,' : 'raterr',
    'Instrumental VEGAMAG magnitude,'    : 'vega',
    'Transformed UBVRI magnitude,'       : 'trans',
    'Magnitude uncertainty,'             : 'err',
    'Chi,'                               : 'chi',
    'Signal-to-noise,'                   : 'snr',
    'Sharpness,'                         : 'sharp',
    'Roundness,'                         : 'round',
    'Crowding,'                          : 'crowd',
    'Photometry quality flag,'           : 'flag',
}



#try:
#    import seaborn as sns; sns.set(style='white', font_scale=1.3)
#except ImportError:
#    print('install seaborn you monster')

# need better way to do this
def name_columns(colfile):
    """Construct a table of column names for dolphot output, with indices
    corresponding to the column number in dolphot output file.

    Inputs
    ------
    colfile : path
        path to file containing dolphot column descriptions

    Returns
    -------
    df : DataFrame
        A table of column descriptions and their corresponding names.
    filters : list
        List of filters included in output
    """
    df = pd.DataFrame(data=np.genfromtxt(colfile, delimiter='. ', dtype=str),
                          columns=['index','desc']).drop('index', axis=1)
    df = df.assign(colnames='')
    # set first 11 column names
    df.loc[:10,'colnames'] = global_columns
    # set rest of column names
    filters_all = []
    for k, v in colname_mappings.items():
        indices = df[df.desc.str.find(k) != -1].index
        desc_split = df.loc[indices,'desc'].str.split(", ")
        # get indices for columns with combined photometry
        indices_total = indices[desc_split.str.len() == 2]
        # get indices for columns with single-frame photometry
        indices_indiv = indices[desc_split.str.len() > 2]
        filters = desc_split.loc[indices_total].str[-1].str.replace("'",'')
        imgnames = desc_split.loc[indices_indiv].str[1].str.split(' ').str[0]
        filters_all.append(filters.values)
        df.loc[indices_total,'colnames'] = filters.str.lower() + '_' + v.lower()
        df.loc[indices_indiv,'colnames'] = imgnames + '_' + v.lower()
    filters_final = np.unique(np.array(filters_all).ravel())
    print('Filters found: {}'.format(filters_final))
    return df, filters_final

def make_cmd(ds, path, targname, red_filter, blue_filter, y_filter, n_err=12,
             #density_kwargs={'f':'log10', 'colormap':'viridis', 'linewidth':2},
             density_kwargs={'f':'log10', 'colormap':'viridis'},
             scatter_kwargs={'c':'k', 'alpha':0.5, 's':1, 'linewidth':2}):
    """Plot a CMD with (blue_filter - red_filter) on the x-axis and
    y_filter on the y-axis.

    Inputs
    ------
    ds : DataFrame
        pandas DataFrame
    red_filter : string
        filter name for "red" filter
    blue_filter : string
        filter name for "blue" filter
    y_filter : string
        filter name for filter on y-axis (usually same as red_filter)
    n_err : int, optional
        number of bins to calculate median photometric errors for
        default: 12
    density_kwargs : dict, optional
        parameters to control density display (e.g. {'f':'log10', 'colormap':'viridis'})
    scatter_kwargs : dict, optional
        parameters to pass to matplotlib.scatter for sparse data

    Returns
    -------
    Nothing

    Outputs
    -------
    some plots dude
    """
    color = '{}-{}'.format(blue_filter.upper(), red_filter.upper())
    blue_vega = '{}_vega'.format(blue_filter)
    red_vega = '{}_vega'.format(red_filter)
    y_vega = '{}_vega'.format(y_filter)
    ylab = '{}'.format(y_filter.upper())
    # ds is a pandas DataFrame here
    ds[color] = ds[blue_vega] - ds[red_vega]
    gst_criteria = ds['{}_gst'.format(red_filter)] & ds['{}_gst'.format(blue_filter)]
    #gst_criteria = ds['({}_gst == True) & ({}_gst == True)'.format(red_filter, blue_filter)]
    name = path + "/" + targname + "_" + blue_filter + "_" + red_filter + "_" + "gst_cmd.png"
    if y_filter not in [blue_filter, red_filter]:
        # idk why you would do this but it's an option
        gst_criteria = gst_criteria & ds['{}_gst'.format(y_filter)]
    # cut dataset down to gst stars (pandas DataFrame)
    ds_gst = ds.loc[gst_criteria].copy()
    # compute plot limits
    xmin = np.nanmin(ds_gst[color].values)
    xmax = np.nanmax(ds_gst[color].values)
    ymin = np.nanmin(ds_gst[y_vega].values)
    ymax = np.nanmax(ds_gst[y_vega].values)
    if xmin < -2.5:
        if "f2" in blue_filter:
            xmin = -3
        else:
            xmin = -2
    if xmax > 10.0:
        if "f2" in blue_filter:
            xmax = 7
        else:
            xmax = 10
    if ymin < 15.0:
        ymin = 15.0
    print(blue_filter, red_filter, " has ", len(ds_gst), " stars in CMD.")

    if len(ds_gst) >= 20000:
        fig, ax = plt.subplots(1, figsize=(7.,5.5))
        plt.rcParams.update({'font.size': 20})
        plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.15)
        # use matplotlib hexbin for dense plots
        data_shape = 200
        # hexbin expects x,y
        hb = ax.hexbin(ds_gst[color].values, ds_gst[y_vega].values,
                       gridsize=data_shape, cmap=density_kwargs.get('colormap','viridis'),
                       bins=None if density_kwargs.get('f', None) is None else 'log')
        fig.colorbar(hb, ax=ax)
        plt.rcParams['axes.linewidth'] = 5
        plt.xticks(np.arange(int(xmin-0.5), int(xmax+0.5), 1.0), fontsize=20)
        plt.yticks(np.arange(int(ymin-0.5), int(ymax+0.5), 1.0), fontsize=20)

        ax.xaxis.set_tick_params(which='minor', direction='in', length=6, width=2, top=True, right=True)
        ax.yaxis.set_tick_params(which='minor', direction='in', length=6, width=2, top=True, right=True)
        ax.xaxis.set_tick_params(direction='in', length=8, width=2, top=True, right=True)
        ax.yaxis.set_tick_params(direction='in', length=8, width=2, top=True, right=True)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(4)
        plt.minorticks_on()
        plt.ylabel(ylab, fontsize=20)
        plt.xlabel(color, fontsize=20)
    else:
        fig, ax = plt.subplots(1, figsize=(7.,5.5), linewidth=5)
        plt.rcParams.update({'font.size': 20})
        plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.15)
        # scatter plot for sparse datasets
        ax.scatter(ds_gst[color].values, ds_gst[y_vega].values, **scatter_kwargs)
        plt.rcParams['axes.linewidth'] = 5
        plt.xticks(np.arange(int(xmin-0.5), int(xmax+0.5), 1.0), fontsize=20)
        plt.yticks(np.arange(int(ymin-0.5), int(ymax+0.5), 1.0), fontsize=20)
        ax.xaxis.set_tick_params(which='minor', direction='in', length=6, width=2, top=True, right=True)
        ax.yaxis.set_tick_params(which='minor', direction='in', length=6, width=2, top=True, right=True)
        ax.xaxis.set_tick_params(direction='in', length=8, width=2, top=True, right=True)
        ax.yaxis.set_tick_params(direction='in', length=8, width=2, top=True, right=True)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(2)
        plt.minorticks_on()
    plt.xlim(int(xmin-0.5), int(xmax+0.5))
    plt.ylim(int(ymin-0.5), int(ymax+0.5))
    plt.ylabel(ylab,fontsize=20)
    plt.xlabel(color,fontsize=20)
    ax.invert_yaxis()
    # compute binned statistics using pandas: bin by y_vega into n_err bins
    bins = np.linspace(ymin, ymax, n_err + 1)
    bin_indices = np.digitize(ds_gst[y_vega].values, bins) - 1
    y_binned = []
    xerr = []
    yerr = []
    x_binned = [xmax * 0.9] * n_err
    for b in range(n_err):
        mask = bin_indices == b
        if np.any(mask):
            y_binned.append(np.nanmean(ds_gst.loc[mask, y_vega].values))
            # x error is sqrt(blue_err^2 + red_err^2)
            xerr.append(np.nanmedian(np.sqrt(ds_gst.loc[mask, '{}_err'.format(blue_filter)].values**2 +
                                             ds_gst.loc[mask, '{}_err'.format(red_filter)].values**2)))
            yerr.append(np.nanmedian(ds_gst.loc[mask, '{}_err'.format(y_filter)].values))
        else:
            y_binned.append(np.nan)
            xerr.append(np.nan)
            yerr.append(np.nan)
    ax.errorbar(x_binned, y_binned, yerr=yerr, xerr=xerr,
                fmt=',', color='k', lw=1.5)
    fig.savefig(name)
    new_dp = wp.DataProduct(my_config, filename=name, 
                             group="proc", data_type="CMD file", subtype="CMD") 


if __name__ == "__main__":
    my_pipe = wp.Pipeline()
    my_job = wp.Job()
    my_job.logprint("processing phot file...")
    my_config = my_job.config
    my_target = my_job.target
    this_event = my_job.firing_event
    my_job.logprint(this_event)
    my_job.logprint(this_event.options)
    my_config = my_job.config
    logpath = my_config.logpath
    procpath = my_config.procpath
    this_dp_id = this_event.options["dp_id"]
    this_dp = wp.DataProduct(int(this_dp_id), group="proc")
    my_job.logprint(
        f"Data Product: {this_dp.filename}\n, Path: {this_dp.target.datapath}\n This DP options{this_dp.options}\n")

    photfile = my_config.procpath  +  "/" + this_dp.filename

    import pandas as pd
    df = pd.read_hdf(photfile, key='data')
    # use pandas DataFrame directly
    ds = df
    #filters = my_config.parameters["det_filters"].split(',')
    colfile = my_config.parameters["colfile"]
    my_job.logprint('Columns file: {}'.format(colfile))
    columns_df, filters = name_columns(colfile)

    waves = []
    for filt in filters:
        pre, suf = filt.split('_')
        wave = suf[1:4]
        if "NIRCAM" in pre:
            wave = 10*int(wave)
        if "WFC3" in pre and "F1" in suf:
            wave = 10*int(wave)
        my_job.logprint(waves)  
        waves.append(int(wave))
    my_job.logprint(waves)  
    sort_inds = np.argsort(waves)
        
    num_all_filters = len(filters)
    for i in range(num_all_filters-1):
       for j in range(num_all_filters-i-1):
           ind2=i+1+j
           my_job.logprint(filters[sort_inds[i]])  
           my_job.logprint(filters[sort_inds[ind2]])  
           try:
               make_cmd(ds, procpath, my_target.name, filters[sort_inds[ind2]].lower(),filters[sort_inds[i]].lower(),filters[sort_inds[ind2]].lower())
           except:
               my_job.logprint(f"{filters[sort_inds[i]]} and {filters[sort_inds[ind2]]} failed")
               continue
       
    next_event = my_job.child_event(
    name="cmds_ready",
    options={"target_id": my_target.target_id}
    )  # next event
    #next_event.fire() 
    #time.sleep(150)
 

    
