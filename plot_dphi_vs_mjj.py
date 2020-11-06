#!/usr/bin/env python

import os
import sys
import re
import matplotlib.colors as colors
from coffea import hist
from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive

pjoin = os.path.join

REBIN = {
    'mjj' : hist.Bin('mjj', r'$M_{jj} \ (GeV)$', 10, 0, 5000),
    'dphi' : hist.Bin('dphi', r'$\Delta \phi (jet,MET)$', 35, 0, 3.5)
}

def plot_mjj(acc,outtag):
    '''Plot QCD MC in signal region as a function of mjj.'''
    variable = 'mjj'
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    mjj_analysis_binning = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])

    h = h.rebin('mjj',mjj_analysis_binning)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        h_temp = h.integrate('dataset', f'QCD_HT_{year}').integrate('region', 'sr_vbf')
        fig, ax = plt.subplots()
        hist.plot1d(h_temp,ax=ax)

        ax.set_yscale('log')
        ax.set_ylim(1e-3,1e6)
        ax.set_title(f'QCD MC SR: {year}')
        ax.get_legend().remove()

        # Save figure
        outpath = pjoin(outdir, f'mjj_qcd_{year}.pdf')
        fig.savefig(outpath)

        plt.close(fig)

def plot_dphi_vs_mjj(acc,outtag):
    '''Plot dphi-mjj 2D distribution for QCD MC, going as input to the QCD estimation.'''
    variable = 'mjj_vs_dphi_qcd'
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    for variable in ['mjj', 'dphi']:
        h = h.rebin(variable, REBIN[variable])

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        # Get QCD MC
        h_temp = h.integrate('dataset', f'QCD_HT_{year}').integrate('region', 'sr_vbf_qcd')
    
        fig, ax = plt.subplots()
        patch_opts = {'norm': colors.LogNorm(vmin=1e-6, vmax=1e5)}
        hist.plot2d(h_temp, xaxis='dphi', ax=ax, patch_opts=patch_opts)
    
        # Save figure
        outpath = pjoin(outdir, f'dphi_vs_mjj_qcd_{year}.pdf')
        fig.savefig(outpath)

        plt.close(fig)

def main():
    inpath = './input/merged_2020-11-06_vbfhinv_03Sep20v7_qcd_estimation_qcd_cr'
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/', '')
    plot_dphi_vs_mjj(acc, outtag)
    plot_mjj(acc, outtag)

if __name__ == '__main__':
    main()