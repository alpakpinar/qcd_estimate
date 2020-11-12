#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

data_err_opts = {
    'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'color':'k'
}

def plot_tf_from_mc(acc, outdir, variable='mjj', ht_startsfrom_200=False):
    '''Calculate and plot the TF from QCD MC.'''
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    if not ht_startsfrom_200:
        h = merge_datasets(h)

    if variable == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200,400,600,900,1200,1500,2000,2750,3500,5000])
        h = h.rebin('mjj', mjj_ax)

    for year in [2017, 2018]:
        if not ht_startsfrom_200:
            _h = h.integrate('dataset', f'QCD_HT_{year}')
        else:
            _h = h[re.compile(f'QCD_HT(?!(50to|100to)).*{year}')].integrate('dataset')
        hnum = _h.integrate('region', 'sr_vbf_qcd_regionB')
        hden = _h.integrate('region', 'sr_vbf_qcd_regionA')

        fig, ax, rax = fig_ratio()
        hist.plot1d(hnum, ax=ax)
        hist.plot1d(hden, ax=ax, clear=False)

        labels = ['Region B', 'Region A']
        ax.legend(labels=labels)
        ax.set_title(f'QCD MC: {year}')
        ax.set_yscale('log')
        ax.set_ylim(1e-3,1e6)

        hist.plotratio(hnum, hden, ax=rax, unc='num', error_opts=data_err_opts)

        rax.grid(True)
        rax.set_ylim(0,1)
        rax.set_ylabel('MC based TF')

        if ht_startsfrom_200:
            outname = f'qcd_mc_regionab_{year}_ht200.pdf'
        else:
            outname = f'qcd_mc_regionab_{year}.pdf'

        outpath = pjoin(outdir, outname)
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc  = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0]

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plot_tf_from_mc(acc, outdir, ht_startsfrom_200=False)
    plot_tf_from_mc(acc, outdir, ht_startsfrom_200=True)

if __name__ == '__main__':
    main()