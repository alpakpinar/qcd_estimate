#!/usr/bin/env python

import os
import sys
import re
import argparse
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint
from stack_plot_qcd_cr import modify_handles_labels, get_fixed_xlabel

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='Path containing merged coffea files.')
    parser.add_argument('--region', help='The region to look at, default is sr_vbf_qcd_cr.', default='sr_vbf_qcd_cr')
    parser.add_argument('--distribution', help='Regex to specify which distribution to plot.', default='.*')
    args = parser.parse_args()
    return args

def plot_jet_fractions_for_eta_slice(acc, distribution, outtag, etaslice, region='sr_vbf_qcd_cr'):
    '''Given the 2D eta/fraction distribution, plot the jet energy fraction for a given eta slice.'''
    # 2D jet eta/NEF distribution for all jets
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Get the relevant eta slice and region
    h = h.integrate('jeteta', etaslice).integrate('region', region)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    for year in [2017, 2018]:
        # Get data and MC 
        h_data = h[f'MET_{year}']
        h_mc = h[ re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}') ]

        fig, ax = plt.subplots()
        hist.plot1d(h_data, ax=ax, overlay='dataset', error_opts=data_err_opts)
        hist.plot1d(h_mc, ax=ax, overlay='dataset', stack=True, clear=False)

        handles, labels = ax.get_legend_handles_labels()
        handles, new_labels = modify_handles_labels(handles, labels)

        ax.legend(handles=handles, labels=new_labels, prop={'size': 10.}, ncol=2)

        ax.set_yscale('log')
        ax.set_ylim(1e-1, 1e4)

        ax.set_title(f'MTR {year}: QCD CR')
        ax.yaxis.set_ticks_position('both')

        # Fix x-label if necessary
        new_xlabel = get_fixed_xlabel(variable='ak4_nef')
        ax.set_xlabel(new_xlabel)

        outdir = f'./output/{outtag}/stack_plot/{region}/em_fracs_by_eta'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        etastarttag = str(etaslice.start).replace('.', '_')
        etaendtag = str(etaslice.stop).replace('.', '_')

        outpath = pjoin(outdir, f'{distribution}_{etastarttag}_{etaendtag}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    args = parse_cli()
    inpath = args.inpath
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/', '')

    etaslices = [
        slice(0,2.5),
        slice(2.5,3.0),
        slice(3.0,5.0)
    ]

    distribtuions = [
        'ak4_eta0_nef0',
        'ak4_eta1_nef1',
        'ak4_eta_nef',
    ]

    for distribtuion in distribtuions:
        if not re.match(args.distribution, distribtuion):
            continue
        for etaslice in etaslices:
            plot_jet_fractions_for_eta_slice(acc, 
                distribution=distribtuion,
                outtag=outtag, 
                etaslice=etaslice, 
                region=args.region
                )

if __name__ == '__main__':
    main()