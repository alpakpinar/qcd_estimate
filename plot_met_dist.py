#!/usr/bin/env python

import os
import sys
import re
import argparse
import numpy as np
from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='Path to merged coffea files.')
    parser.add_argument('--region', help='Region to look at, default is inclusive.', default='inclusive')
    args = parser.parse_args()
    return args

def get_title(region, year):
    if region == 'inclusive':
        return f'QCD MC Inclusive: {year}'
    return f'QCD MC {region}: {year}'

def plot_met_dist(acc, outtag, region='inclusive', split_by_ht=False):
    '''Plot MET distribution for QCD.'''
    variable = 'met'
    acc.load(variable)
    h = acc[variable]

    # Output directory to save plots
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    
    # Rebin MET
    met_ax = hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,20)))
    h = h.rebin('met', met_ax)

    if not split_by_ht:
        h = merge_datasets(h)

    h = h.integrate('region', region)
    for year in [2017, 2018]:
        fig, ax = plt.subplots()
        if not split_by_ht:
            _h = h.integrate('dataset', f'QCD_HT_{year}')
            hist.plot1d(_h, ax=ax)
            ax.get_legend().remove()
        else:
            _h = h[re.compile(f'QCD_HT.*{year}')]
            hist.plot1d(_h, ax=ax, overlay='dataset')

            # Handle legend
            new_labels = []
            _, labels = ax.get_legend_handles_labels()
            for label in labels:
                new_label = label.replace('QCD_','').replace(f'_{year}', '')
                new_label = re.sub('-(mg|MLM)', '', new_label)
                new_labels.append(new_label)

            ax.legend(labels=new_labels)

        ax.set_title( get_title(region, year) )

        suffix = '_htsplit' if split_by_ht else ''

        outpath = pjoin(outdir, f'met_dist_{region}_{year}{suffix}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    args = parse_cli()
    inpath = args.inpath
    region = args.region
    
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag  = re.findall('merged_.*', inpath)[0]

    plot_met_dist(acc, outtag, region=args.region, split_by_ht=False)
    plot_met_dist(acc, outtag, region=args.region, split_by_ht=True)

if __name__ == '__main__':
    main()