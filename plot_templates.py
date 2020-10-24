#!/usr/bin/env python

import os
import sys
import re
import uproot
import mplhep as hep
import numpy as np
from matplotlib import pyplot as plt

pjoin = os.path.join

def plot_qcd(inpath, fit='nominal', binning='nom', region='sr_vbf_qcd', cr_only=False):
    '''Plot QCD templates in SR and QCD CR.'''
    fpath = pjoin(inpath, f'templates_{region}_{fit}_bin_{binning}.root')

    # Output directory to save plots
    outdir = pjoin(inpath, 'templates')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    infile = uproot.open(fpath)
    for year in [2017, 2018]:
        # Plot QCD in SR and CR 
        h_sr = infile[f'{region}_{year}_sr_qcd']
        h_cr = infile[f'{region}_{year}_cr_qcd']

        fig, ax = plt.subplots()
        hep.histplot(h_cr.values, h_cr.edges, yerr=np.sqrt(h_cr.variances), label='CR')
        if not cr_only:
            hep.histplot(h_sr.values, h_sr.edges, yerr=np.sqrt(h_sr.variances), label='SR')

        ax.set_xlabel(r'$M_{jj} \ (GeV)$')
        ax.set_ylabel('Counts')
        ax.set_title(f'QCD MC: {year}')
        ax.legend()

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1e8)

        # Save figure
        if cr_only:
            outname = f'qcd_template_cr_bin_{binning}_{year}.pdf'
        else:
            outname = f'qcd_templates_bin_{binning}_{year}.pdf'
        outpath = pjoin(outdir, outname)
        fig.savefig(outpath)
        print(f'MSG% File saved: {outpath}')

def main():
    for binning in ['nom']:
        for region in ['sr_vbf_qcd_recoil_200', 'sr_vbf_qcd_recoil_230']:
            # Input path for the template root files
            inpath = f'output/merged_2020-10-23_vbfhinv_03Sep20v7_qcd_estimation_loose_recoil_regions/{region}'
        
            try:
                plot_qcd(inpath, fit='nominal', binning=binning, region=region)
                plot_qcd(inpath, fit='nominal', binning=binning, cr_only=True, region=region)
            except KeyError:
                print(f'Could not find binning: {binning}, skipping')
                continue
    
if __name__ == '__main__':
    main()