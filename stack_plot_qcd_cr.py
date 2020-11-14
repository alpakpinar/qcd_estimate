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

legend_labels = {
    'DY.*' : "QCD Z$\\rightarrow\\ell\\ell$",
    'EWKZ.*ZToLL.*' : "EWK Z$\\rightarrow\\ell\\ell$",
    'WN*J.*LNu.*' : "QCD W$\\rightarrow\\ell\\nu$",
    'EWKW.*LNu.*' : "EWK W$\\rightarrow\\ell\\nu$",
    'ZJetsToNuNu.*.*' : "QCD Z$\\rightarrow\\nu\\nu$",
    'EWKZ.*ZToNuNu.*' : "EWK Z$\\rightarrow\\nu\\nu$",
    'DY.*' : "QCD Z$\\rightarrow\\ell\\ell$",
    'WN*J.*LNu.*' : "QCD W$\\rightarrow\\ell\\nu$",
    'ZJetsToNuNu.*.*' : "QCD Z$\\rightarrow\\nu\\nu$",
    'QCD.*' : "QCD",
    'Top.*' : "Top quark",
    'Diboson.*' : "WW/WZ/ZZ",
    'MET|Single(Electron|Photon|Muon)|EGamma.*' : "Data"
}

binnings = {
    'ak4_pt0' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,1000,20)) ),
    'ak4_pt1' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,800,20)) ),
    'mjj' : hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
}
            
def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='Path containing merged coffea files.')
    parser.add_argument('--region', help='The region to be plotted, default is QCD CR.', default='sr_vbf_qcd_cr')
    parser.add_argument('--distribution', help='Regex for the distributions to be plotted.', default='.*')
    args = parser.parse_args()
    return args

def modify_handles_labels(handles, labels):
    '''Helper function to update legend labels and aesthetics.'''
    new_labels = []
    for handle, label in zip(handles, labels):
        if not 'MET' in label:
            handle.set_linestyle('-')
            handle.set_edgecolor('k')

        for k, v in legend_labels.items():
            if re.match(k, label):
                l = v
        new_labels.append(l if l else label)

    return handles, new_labels

def fix_xlabel(ax, variable):
    xlabels_to_fix = {
        'ak4_nef.*' : 'Neutral EM Energy Fraction',
        'ak4_nhf.*' : 'Neutral Hadron Energy Fraction',
        'ak4_chf.*' : 'Charged Hadron Energy Fraction',
        'ak4_eta0' : r'Leading Jet $\eta$',
        'ak4_eta1' : r'Trailing Jet $\eta$',
        'dphitkpf' : r'$\Delta\phi(TkMET, PFMET)$',
    }
    for key, xlabel in xlabels_to_fix.items():
        if re.match(key, variable):
            ax.set_xlabel(xlabel)

    return ax

def stack_plot_qcd_cr(acc, outtag, variable='detajj', region='sr_vbf_qcd_cr'):
    '''Create a stack plot for QCD CR.'''
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Get the QCD CR
    h = h.integrate('region', region)

    # Rebin, if necessary
    try:
        newax = binnings[variable]
        h = h.rebin(h.axis(newax.name), newax)
    except KeyError:
        pass

    for year in [2017, 2018]:
        h_data = h[f'MET_{year}']
        h_mc = h[re.compile(f'(ZJetsToNuNu.*|EW.*|QCD.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}')]
        
        data_err_opts = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
        }

        if variable == 'vecdphi':
            overflow='over'
        else:
            overflow='none'

        fig, ax = plt.subplots()
        hist.plot1d(h_data, ax=ax, overlay='dataset', error_opts=data_err_opts, overflow=overflow)
        hist.plot1d(h_mc, ax=ax, stack=True, overlay='dataset', clear=False, overflow=overflow)

        handles, labels = ax.get_legend_handles_labels()
        handles, new_labels = modify_handles_labels(handles, labels)

        ax.legend(handles=handles, labels=new_labels, prop={'size': 10.}, ncol=2)

        ax.set_yscale('log')
        ax.set_ylim(1e-1, 1e4)

        ax.set_title(f'MTR {year}: QCD CR')
        ax.yaxis.set_ticks_position('both')

        if variable in ['ak4_eta0', 'ak4_eta1']:
            ax.axvline(x=3.0, ymin=0, ymax=1, color='red', lw=2)
            ax.axvline(x=-3.0, ymin=0, ymax=1, color='red', lw=2)

        # Fix x-label if necessary
        ax = fix_xlabel(ax, variable)

        outdir = f'./output/{outtag}/stack_plot/{region}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        outpath = pjoin(outdir, f'stack_plot_{variable}_{year}.pdf')
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

    variables = [
        'detajj',
        'dphijj',
        'mjj',
        'ak4_pt0',
        'ak4_eta0',
        'ak4_nef0',
        'ak4_nhf0',
        'ak4_chf0',
        'ak4_pt1',
        'ak4_eta1',
        'ak4_nef1',
        'ak4_nhf1',
        'ak4_chf1',
        'vecb',
        'vecdphi',
        'dphitkpf'
    ]

    for variable in variables:
        if not re.match(args.distribution, variable):
            continue
        try:
            stack_plot_qcd_cr(acc, outtag, region=args.region, variable=variable)
        except KeyError:
            print(f'Distribution not found: {variable}, skipping')

if __name__ == '__main__':
    main()
