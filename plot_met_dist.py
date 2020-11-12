#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive

pjoin = os.path.join

def plot_met_dist(acc, outtag, region='inclusive'):
    '''Plot MET distribution for QCD.'''
    variable = 'met'
    acc.load(variable)
    h = acc[variable]

    # Output directory to save plots
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Rebin MET
    met_ax = hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,20)))
    h = h.rebin('met', met_ax)

    h = h.integrate('region', region)
    for year in [2017, 2018]:
        _h = h.integrate('dataset', f'QCD_HT_{year}')
        fig, ax = plt.subplots()
        hist.plot1d(_h, ax=ax)

        ax.get_legend().remove()
        ax.set_title(f'QCD MC Inclusive: {year}')

        outpath = pjoin(outdir, f'met_dist_{region}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag  = re.findall('merged_.*', inpath)[0]

    plot_met_dist(acc, outtag)

if __name__ == '__main__':
    main()