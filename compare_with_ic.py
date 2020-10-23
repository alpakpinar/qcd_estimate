#!/usr/bin/env python

import os
import sys
import re
import uproot
import mplhep as hep
import numpy as np
import matplotlib.ticker
from matplotlib import pyplot as plt
from bucoffea.plot.util import fig_ratio

pjoin = os.path.join

def compare_bu_ic_qcd_estimation(ic_inputdir, jobtag):
    '''Make comparison plots of QCD estimations of BU and IC.'''
    bu_inputfile = './qcdestimate_sr_vbf_qcd.root'
    f_bu = uproot.open(bu_inputfile)

    outdir = f'./output/{jobtag}/bu_ic_comp'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        f_ic = uproot.open(pjoin(ic_inputdir, f'vbf_qcd_nckw_ws_{year}.root'))[f'category_vbf_{year}']
        h_ic = f_ic['rebin_QCD_hist_counts']
        
        h_bu = f_bu[f'qcd_vbf_{year}']

        fig, ax, rax = fig_ratio()
        hep.histplot(h_ic.values, h_ic.edges, ax=ax, binwnorm=True, label='IC')
        hep.histplot(h_bu.values, h_bu.edges, ax=ax, binwnorm=True, label='BU')

        # Compare with QCD MC
        f_qcdmc = uproot.open(f'./output/{jobtag}/templates_sr_vbf_qcd_nominal_bin_nom.root')
        h_qcdmc = f_qcdmc[f'sr_vbf_qcd_{year}_sr_qcd']

        hep.histplot(h_qcdmc.values, h_qcdmc.edges, yerr=np.sqrt(h_qcdmc.variances), ax=ax, label='QCD MC in SR', histtype='errorbar', binwnorm=True, color='k')

        ax.set_yscale('log')
        ax.set_ylim(1e-5, 1e4)
        ax.set_ylabel('Events / GeV')
        ax.set_title(f'MTR {year}')
        ax.legend(title='QCD Templates')

        loc = matplotlib.ticker.LogLocator(base=10, numticks=15)
        ax.yaxis.set_major_locator(loc)

        # Plot ratio
        centers = 0.5 * np.sum(h_ic.bins, axis=1)

        rax.plot(centers, h_bu.values / h_ic.values, marker='o', ls='', color='k')
        rax.set_xlabel(r'$M_{jj} \ (GeV)$')
        rax.set_ylabel('BU / IC')
        rax.set_ylim(0,2)
        rax.grid(True)

        # Save figure
        outpath = pjoin(outdir, f'bu_ic_comp_qcd_{year}.pdf')
        fig.savefig(outpath)
        print(f'File saved: {outpath}')

def main():
    # Path to input files from IC
    ic_inputdir = './input/from_ic'

    jobtag = 'merged_2020-10-22_vbfhinv_03Sep20v7_qcd_estimation'

    compare_bu_ic_qcd_estimation(ic_inputdir, jobtag)

if __name__ == '__main__':
    main()