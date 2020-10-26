#!/usr/bin/env python

import os
import sys
import re
import uproot
import mplhep as hep
import numpy as np
from matplotlib import pyplot as plt
from data_driven_qcd import histdiff
from klepto.archives import dir_archive
from coffea import hist
from pprint import pprint

pjoin = os.path.join

def calculate_integral(h):
    '''Calculate the integral of the mjj distribution.'''
    bin_widths = np.diff(h.edges)
    bin_vals = h.values
    return np.sum(bin_vals * bin_widths)

def plot_qcd(inpath, fit='nominal', binning='nom', cr_only=False):
    '''Plot QCD MC templates in SR and QCD CR.'''
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

        print(f'Year: {year}')
        integral_sr = calculate_integral(h_sr)
        print(f'Integral of QCD MC over SR: {integral_sr}')
        integral_cr = calculate_integral(h_cr)
        print(f'Integral of QCD MC over CR: {integral_cr}')

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
        plt.close(fig)
        print(f'MSG% File saved: {outpath}')

def plot_data_driven_qcd_in_cr(inpath, fit='nominal', binning='nom'):
    '''In the QCD CR, plot the data driven QCD estimation, e.g. Data - Non-QCD bkg'''
    fpath = pjoin(inpath, f'templates_sr_vbf_qcd_{fit}_bin_{binning}.root')

    # Output directory to save plots
    outdir = pjoin(inpath, 'templates')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    infile = uproot.open(fpath)
    for year in [2017, 2018]:
        # Read the data and non-QCD histograms from file
        h_cr_data = infile[f'sr_vbf_qcd_{year}_cr_data']
        h_cr_nonqcd = infile[f'sr_vbf_qcd_{year}_cr_nonqcd']

        # Take the difference of the two to get the data-driven QCD estimate in CR
        qcd_est_sumw, qcd_est_sumw2 = histdiff(h_cr_data, h_cr_nonqcd)

        # Plot these
        fig, ax = plt.subplots()
        hep.histplot(h_cr_data.values, h_cr_data.edges, yerr=np.sqrt(h_cr_data.variances), label='Data')
        hep.histplot(h_cr_nonqcd.values, h_cr_nonqcd.edges, yerr=np.sqrt(h_cr_nonqcd.variances), label='Non-QCD bkg')

        # Plot the difference 
        hep.histplot(qcd_est_sumw, h_cr_data.edges, yerr=np.sqrt(qcd_est_sumw2), histtype='errorbar', color='k', label='QCD Estimate in CR')

        ax.set_xlabel(r'$M_{jj} \ (GeV)$')
        ax.set_ylabel('Counts')
        ax.legend()
        ax.set_title(f'QCD Prediction in CR: {year}')

        ax.set_yscale('log')
        ax.set_ylim(1e-2, 1e8)

        # Save figure
        outpath = pjoin(outdir, f'qcd_in_cr_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'MSG% File saved: {outpath}')

def plot_templates_split_by_ht(inpath, region='sr'):
    '''Plot QCD templates in SR or CR, split by the individual HT bins.'''
    acc = dir_archive(inpath)
    acc.load('sumw')
    variable = 'mjj_vs_dphi_qcd'
    acc.load(variable)
    h = acc[variable]

    from bucoffea.plot.util import merge_extensions, scale_xs_lumi
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)

    # Get the relevant dphi values
    if region == 'sr':
        dphi_slice = slice(0.5, None)
    else:
        dphi_slice = slice(0., 0.5)
    
    h = h.integrate('dphi', dphi_slice).integrate('region', 'sr_vbf_qcd')
    # Rebin mjj
    mjj_bin = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200,400,600,900,1200,1500,2000,2750,3500,5000])
    h = h.rebin('mjj', mjj_bin)

    # Output directory for plots
    outdir = pjoin( inpath.replace('input', 'output'), 'split_by_ht' )
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        # Get QCD MC datasets for the given year
        _h = h[re.compile(f'QCD_HT.*{year}')]

        # Plot for each HT bin
        fig, ax = plt.subplots()
        hist.plot1d(_h, ax=ax, overlay='dataset')

        ax.set_yscale('log')
        ax.set_ylim(1e-4, 1e8)
        if region == 'sr':
            ax.set_title(r'$\Delta\phi(jet,MET) > 0.5$')
        else:
            ax.set_title(r'$\Delta\phi(jet,MET) \leq 0.5$')

        # Make the legend labels shorter
        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            new_label = label.split('_')[1].split('-')[0]
            handle.set_label(new_label)

        ax.legend(handles=handles, prop={'size' : 10.})

        # Save figure
        outpath = pjoin(outdir, f'qcd_split_by_ht_{region}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

def main():
    # Input path for the template root files
    inpath = 'output/merged_2020-10-22_vbfhinv_03Sep20v7_qcd_estimation/'

    for binning in ['nom', 'alt1', 'alt2', 'alt3', 'alt4']:
        try:
            plot_qcd(inpath, fit='nominal', binning=binning)
            plot_qcd(inpath, fit='nominal', binning=binning, cr_only=True)

            plot_data_driven_qcd_in_cr(inpath, fit='nominal', binning=binning)
        except KeyError:
            print(f'Could not find binning: {binning}, skipping')
            continue
    
    inpath_for_klepto = 'input/merged_2020-10-22_vbfhinv_03Sep20v7_qcd_estimation/'
    plot_templates_split_by_ht(inpath_for_klepto, region='sr')
    plot_templates_split_by_ht(inpath_for_klepto, region='cr')

if __name__ == '__main__':
    main()