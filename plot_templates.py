#!/usr/bin/env python

import os
import sys
import re
import uproot
import mplhep as hep
import numpy as np
from matplotlib import pyplot as plt

pjoin = os.path.join

def plot_qcd(inpath, fit='nominal', binning='nom'):
    '''Plot QCD templates in SR and QCD CR.'''
    fpath = pjoin(inpath, f'templates_sr_vbf_qcd_{fit}_bin_{binning}.root')

    # Output directory to save plots
    outdir = pjoin(inpath, 'templates')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    infile = uproot.open(fpath)
    for year in [2017, 2018]:
        # Plot QCD in SR and CR 
        h_sr = infile[f'sr_vbf_qcd_{year}_sr_qcd']
        h_cr = infile[f'sr_vbf_qcd_{year}_cr_qcd']

        fig, ax = plt.subplots()
        hep.histplot(h_sr.values, h_sr.edges, yerr=np.sqrt(h_sr.variances), label='SR')
        hep.histplot(h_cr.values, h_cr.edges, yerr=np.sqrt(h_cr.variances), label='CR')

        ax.set_xlabel(r'$M_{jj} \ (GeV)$')
        ax.set_ylabel('Counts')
        ax.set_title(f'QCD MC: {year}')
        ax.legend()

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1e8)

        # Save figure
        outpath = pjoin(outdir, f'qcd_templates_{year}.pdf')
        fig.savefig(outpath)

def main():
    # Input path for the template root files
    inpath = 'output/merged_2020-10-22_vbfhinv_03Sep20v7_qcd_estimation/'

    plot_qcd(inpath, fit='nominal', binning='nom')

if __name__ == '__main__':
    main()