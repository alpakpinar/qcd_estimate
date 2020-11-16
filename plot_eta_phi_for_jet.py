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

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='Path containing merged coffea files.')
    parser.add_argument('--region', help='The region to look at, default is sr_vbf_qcd_cr_detajj.', default='sr_vbf_qcd_cr_detajj')
    args = parser.parse_args()
    return args

def plot_eta_phi_map(acc, outtag, distribution, region):
    '''Plot 2D eta/phi map for leading or trailing jet.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Integrate the region
    h = h.integrate('region', region)

    # Output directory to save plots
    outdir = f'./output/{outtag}/jet_eta_phi/{region}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    for year in [2017, 2018]:
        # Get data and MC 
        h_data = h.integrate('dataset', f'MET_{year}')
        h_mc = h.integrate('dataset', re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}'))

        # Plot 2D map for data / MC!
        fig, ax = plt.subplots()
        etaedges = h_data.axis('jeteta').edges()
        phiedges = h_data.axis('jetphi').edges()

        ratio = h_data.values()[()] / h_mc.values()[()]

        pc = ax.pcolormesh(etaedges, phiedges, ratio.T)
        fig.colorbar(pc, ax=ax, label='Data / MC')

        ax.set_xlabel(r'Jet $\eta$')
        ax.set_ylabel(r'Jet $\phi$')

        # Save figure
        outpath = pjoin(outdir, f'{distribution}_{year}.pdf')
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

    plot_eta_phi_map(acc, outtag,
        distribution='ak4_eta0_phi0',
        region=args.region
    )

if __name__ == '__main__':
    main()