#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from coffea import hist
from matplotlib import pyplot as plt
from klepto.archives import dir_archive

pjoin = os.path.join

def stack_plot_qcd_cr(acc, outtag, variable='detajj'):
    '''Create a stack plot for QCD CR.'''
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Get the QCD control region
    h = h.integrate('region', 'sr_vbf_qcd_cr')

    for year in [2017, 2018]:
        # Get data and MC
        h_data = h[f'MET_{year}']
        h_mc = h[re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}')]

        fig, ax = plt.subplots()
        hist.plot1d(h_data, ax=ax, overlay='dataset')
        hist.plot1d(h_mc, ax=ax, overlay='dataset', clear=False)

        # Save figure
        outdir = f'./output/{outtag}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'stack_plot_qcd_cr_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    # UPDATE!
    inpath = './input/merged_2020-10-27_vbfhinv_03Sep20v7_qcd_estimation_very_loose_recoil_regions_detajj_cat'
    outtag = inpath.split('/')[-1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    stack_plot_qcd_cr(acc, outtag)

if __name__ == '__main__':
    main()