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
    _temp = task.mask(source="*", name="fake_hdf5_ready", value="*")




#try:
#    import seaborn as sns; sns.set(style='white', font_scale=1.3)
#except ImportError:
#    print('install seaborn you monster')

# need better way to do this

def name_columns(colfile):
    df = pd.DataFrame(data=np.genfromtxt(colfile, delimiter='. ', dtype=str),
                          columns=['index','desc']).drop('index', axis=1)
    df = df.assign(colnames='')
    # set first 11 column names
    df.loc[:10,'colnames'] = ['ext','chip','x','y','chi_gl','snr_gl','sharp_gl',
                              'round_gl','majax_gl','crowd_gl','objtype_gl']
    filters_all = []
    colname_mappings = {
        'counts,': 'count', 'sky level,': 'sky', 'Normalized count rate,': 'rate',
        'Normalized count rate uncertainty,': 'raterr', 'Instrumental VEGAMAG magnitude,': 'vega',
        'Transformed UBVRI magnitude,': 'trans', 'Magnitude uncertainty,': 'err', 'Chi,': 'chi',
        'Signal-to-noise,': 'snr', 'Sharpness,': 'sharp', 'Roundness,': 'round',
        'Crowding,': 'crowd', 'Photometry quality flag,': 'flag'
    }
    for k, v in colname_mappings.items():
        indices = df[df.desc.str.find(k) != -1].index
        desc_split = df.loc[indices,'desc'].str.split(', ')
        indices_total = indices[desc_split.str.len() == 2]
        indices_indiv = indices[desc_split.str.len() > 2]
        filters = desc_split.loc[indices_total].str[-1].str.replace("'", '')
        imgnames = desc_split.loc[indices_indiv].str[1].str.split(' ').str[0]
        filters_all.append(filters.values)
        df.loc[indices_total,'colnames'] = filters.str.lower() + '_' + v.lower()
        df.loc[indices_indiv,'colnames'] = imgnames + '_' + v.lower()
    filters_final = np.unique(np.array(filters_all).ravel())
    print('Filters found: {}'.format(filters_final))
    return df, filters_final

def make_resid_plot(my_job, ds, path, targname, filter, n_err=12,
                    density_kwargs={'f':'log10', 'colormap':'viridis'},
                    scatter_kwargs={'c':'k', 'alpha':0.5, 's':1, 'linewidth':2}):
    """Plot residuals (Out - In) vs input magnitude for a filter.

    ds is expected to be a pandas DataFrame.
    """
    vega = '{}_vega'.format(filter)
    # find input mag column (ends with _magin)
    check = 0
    incolname = None
    for colname in ds.columns:
        if "magin" not in colname or check > 0:
            continue
        image = colname.split("_magin")[0]
        # prefer using job.config if available
        try:
            cfg_id = my_job.config.id
        except Exception:
            cfg_id = None
        try:
            imdp = wp.DataProduct.select(config_id=cfg_id, filename=image + ".fits", group="proc")
            camera = imdp.options.get("camera", "")
            filt = imdp.options.get("filter", "")
            camfilt = camera + "_" + filt
            if camfilt in filter:
                incolname = colname
                check += 1
        except Exception:
            # fallback: accept the first magin column
            incolname = colname
            check += 1
    if incolname is None:
        raise ValueError(f"No input mag column found for filter {filter}")

    xlab = f"{filter.upper()} IN"
    ylab = "Out - In"
    gst_criteria = ds['{}_gst'.format(filter)]
    name = path + "/" + targname + "_" + filter + "_" + "gst_asts.png"
    ds_gst = ds.loc[gst_criteria].copy()

    xmin = np.nanmin(ds_gst[incolname].values)
    xmax = np.nanmax(ds_gst[incolname].values)
    ymin = -1.0
    ymax = 1.0
    diff = ds_gst[vega].values - ds_gst[incolname].values
    my_job.logprint(f"{filter} has {len(ds_gst)}/{len(ds)} stars recovered.")

    fig, ax = plt.subplots(1, figsize=(7., 5.5))
    plt.rcParams.update({'font.size': 20})
    plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.15)

    data_shape = 200
    if len(ds_gst) >= 50000:
        hb = ax.hexbin(ds_gst[incolname].values, diff,
                       gridsize=data_shape, cmap=density_kwargs.get('colormap', 'viridis'), bins=None)
        fig.colorbar(hb, ax=ax)
    else:
        ax.scatter(ds_gst[incolname].values, diff, **scatter_kwargs)

    plt.xticks(np.arange(int(xmin - 0.5), int(xmax + 0.5), 1.0), fontsize=20)
    plt.yticks(np.arange(int(ymin - 0.5), int(ymax + 0.5), 1.0), fontsize=20)
    ax.xaxis.set_tick_params(which='minor', direction='in', length=6, width=2, top=True, right=True)
    ax.yaxis.set_tick_params(which='minor', direction='in', length=6, width=2, top=True, right=True)
    ax.xaxis.set_tick_params(direction='in', length=8, width=2, top=True, right=True)
    ax.yaxis.set_tick_params(direction='in', length=8, width=2, top=True, right=True)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2)
    plt.minorticks_on()
    plt.xlim(int(xmin - 0.5), int(xmax + 0.5))
    plt.ylim(int(ymin - 0.5), int(ymax + 0.5))
    plt.ylabel(ylab, fontsize=20)
    plt.xlabel(xlab, fontsize=20)
    ax.invert_yaxis()

    # binned statistics along input mag
    bins = np.linspace(xmin, xmax, n_err + 1)
    bin_idx = np.digitize(ds_gst[incolname].values, bins) - 1
    x_binned = []
    y_binned = []
    yerr = []
    for b in range(n_err):
        mask = bin_idx == b
        if np.any(mask):
            x_binned.append(np.nanmean(ds_gst.loc[mask, incolname].values))
            y_med = np.nanmedian(diff[mask])
            y_binned.append(y_med)
            yerr.append(np.nanmedian(np.abs(diff[mask] - y_med)))
        else:
            x_binned.append(np.nan)
            y_binned.append(np.nan)
            yerr.append(np.nan)

    ax.errorbar(x_binned, y_binned, yerr=yerr, fmt=',', color='k', lw=1.5)
    fig.savefig(name)
    new_dp = wp.DataProduct(my_config, filename=name, group="proc", data_type="CMD file", subtype="CMD")
    my_job.logprint("processing phot file...")
    my_config = my_job.config
    my_target = my_job.target
    this_event = my_job.firing_event
    my_job.logprint(this_event)
    my_job.logprint(this_event.options)
    logpath = my_config.logpath
    procpath = my_config.procpath
    this_dp_id = this_event.options["dp_id"]
    this_dp = wp.DataProduct(int(this_dp_id), group="proc")
    my_job.logprint(
        f"Data Product: {this_dp.filename}\n, Path: {this_dp.target.datapath}\n This DP options{this_dp.options}\n")

    photfile = this_dp.filename

    import pandas as pd
    df = pd.read_hdf(photfile, key='data')
    ds = df
    # get column mapping and filters
    colfile = my_config.parameters.get("colfile")
    if colfile is None:
        raise ValueError("colfile parameter not set in config")
    my_job.logprint('Columns file: {}'.format(colfile))
    columns_df, filters = name_columns(colfile)

    waves = []
    for filt in filters:
        pre, suf = filt.split('_')
        wave = suf[1:4]
        if "NIRCAM" in pre:
            wave = 10 * int(wave)
        if "WFC3" in pre and "F1" in suf:
            wave = 10 * int(wave)
        waves.append(int(wave))
    sort_inds = np.argsort(waves)

    num_all_filters = len(filters)
    for i in range(num_all_filters - 1):
        for j in range(num_all_filters - i - 1):
            ind2 = i + 1 + j
            my_job.logprint(filters[sort_inds[i]])
            my_job.logprint(filters[sort_inds[ind2]])
            try:
                make_resid_plot(my_job, ds, procpath, my_target.name, filters[sort_inds[ind2]].lower())
            except Exception as e:
                my_job.logprint(f"An error occurred: {e}")
                my_job.logprint(f"{filters[sort_inds[i]]} and {filters[sort_inds[ind2]]} failed")
                continue

    next_event = my_job.child_event(name="cmds_ready", options={"target_id": my_target.target_id})
    next_event.fire()
    time.sleep(150)
 

    
