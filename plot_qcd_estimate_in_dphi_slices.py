#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from matplotlib import pyplot as plt
from klepto.archives import dir_archive

pjoin = os.path.join

def plot_qcd_estimate_in_dphi_slices(acc, outdir, year, slices=[ slice(0,0.1), slice(0.1,0.4), slice(0.4,0.5) ]):
    '''Plot QCD estimate in CR in different dphi slices, as specified by the slices argument.'''
    variable = 'mjj_vs_dphi_qcd'
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Rebin mjj
    mjj_bin = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500])
    h = h.rebin('mjj', mjj_bin)

    h = h.integrate('region', 'sr_vbf_qcd')

    fig, ax, rax = fig_ratio()
    h_qcd_list = []
    for dphi_slice in slices:
        _h = h.integrate('dphi', dphi_slice)

        # For each dphi slice, calculate the QCD estimation by getting data - non-QCD backgrounds
        h_qcd = _h[f'MET_{year}'].integrate('dataset')
        h_mc = _h[re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')].integrate("dataset")

        h_mc.scale(-1)
        h_qcd.add(h_mc)

        hist.plot1d(h_qcd, ax=ax, clear=False, density=True)
        h_qcd_list.append(h_qcd)

    labels = [
        r'$\Delta\phi_{jm} < 0.1$',
        r'$0.1 < \Delta\phi_{jm} < 0.2$',
        r'$0.2 < \Delta\phi_{jm} < 0.5$',
    ]

    # Fix the labels
    ax.legend(labels=labels)

    ax.set_yscale('log')
    ax.set_ylim(1e-6,1e0)
    ax.set_xlabel('')
    ax.set_ylabel('Normalized Counts')
    ax.set_title(f'MTR {year}: QCD Estimate in CR')

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 12.,
    }

    # Ratios
    base_qcd = h_qcd_list[0]
    for idx, h_qcd in enumerate(h_qcd_list[1:]):
        data_err_opts['color'] = f'C{idx+1}'
        hist.plotratio(
            h_qcd,
            base_qcd,
            ax=rax,
            unc='num',
            error_opts=data_err_opts,
            clear=False
        )

    rax.grid(True)
    rax.set_ylabel(r'Ratio to $\Delta\phi_{jm} < 0.1$')
    rax.set_ylim(0,2)    
    rax.set_xlabel(r'$M_{jj} \ (GeV)$')

    rax.axhline(1, 0, 1, color='k')

    # Save figure
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'qcd_estimate_by_dphi_{year}.pdf')

    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = './input/merged_2020-11-06_vbfhinv_03Sep20v7_qcd_estimation_qcd_cr'
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    # Output directory to save plots
    outdir = inpath.replace('input', 'output')

    dphi_slices = [
        slice(0, 0.1),
        slice(0.1, 0.2),
        slice(0.2, 0.5)
    ]

    for year in [2017, 2018]:
        plot_qcd_estimate_in_dphi_slices(acc, outdir, year=year, slices=dphi_slices)

if __name__ == '__main__':
    main()

