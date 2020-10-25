#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from bucoffea.plot.style import matplotlib_rc
from coffea import hist
from matplotlib import pyplot as plt
from klepto.archives import dir_archive

pjoin = os.path.join

matplotlib_rc()

def plot_recoil(acc, outtag, region='sr_vbf_qcd_recoil_200', year=2017):
    variable = 'recoil'
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Get QCD MC in the requested region
    h = h.integrate('region', region).integrate('dataset', f'QCD_HT_{year}')
    
    # Plot the recoil distributon
    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)
    ax.set_xlim(200,1000)

    ylim = ax.get_ylim()
    ax.plot([250.,250.], ylim, color='red', lw=2)
    ax.set_ylim(ylim)

    ax.get_legend().remove()

    ax.set_title(f'QCD MC: {year}')

    # Save figure
    outdir = f'./output/{outtag}/{region}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'qcd_mc_recoil_dist_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    # Path to merged coffea files
    inpath = './input/merged_2020-10-23_vbfhinv_03Sep20v7_qcd_estimation_loose_recoil_regions'
    acc = dir_archive(inpath)

    acc.load('sumw')
    acc.load('sumw2')

    outtag = inpath.split('/')[-1]

    for year in [2017, 2018]:
        for region in ['sr_vbf_qcd_recoil_200', 'sr_vbf_qcd_recoil_230']:
            plot_recoil(acc, outtag, region=region, year=year)

if __name__ == '__main__':
    main()
