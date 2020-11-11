#!/usr/bin/env python

import os
import sys
import re
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

def stack_plot_qcd_cr(acc, outtag, variable='detajj', region='sr_vbf_qcd_cr'):
    '''Create a stack plot for QCD CR.'''
    acc.load(variable)
    h = acc[variable]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Get the QCD CR
    h = h.integrate('region', region)

    for year in [2017, 2018]:
        h_data = h[f'MET_{year}']
        h_mc = h[re.compile(f'(ZJetsToNuNu.*|EW.*|QCD.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}')]
        
        data_err_opts = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
        }

        fig, ax = plt.subplots()
        hist.plot1d(h_data, ax=ax, overlay='dataset', error_opts=data_err_opts)
        hist.plot1d(h_mc, ax=ax, stack=True, overlay='dataset', clear=False)

        handles, labels = ax.get_legend_handles_labels()
        new_labels = []
        for handle, label in zip(handles, labels):
            if not 'MET' in label:
                handle.set_linestyle('-')
                handle.set_edgecolor('k')

            for k, v in legend_labels.items():
                if re.match(k, label):
                    l = v
            new_labels.append(l if l else label)

        ax.legend(handles=handles, labels=new_labels, prop={'size': 10.})

        ax.set_yscale('log')
        ax.set_ylim(1e-1, 1e4)

        ax.set_title(f'MTR {year}: QCD CR')

        outdir = f'./output/{outtag}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        outpath = pjoin(outdir, f'stack_plot_{variable}_{region}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    inpath = './input/merged_2020-11-11_vbfhinv_03Sep20v7_qcd_estimation_qcd_cr_detajj'
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = inpath.split('/')[-1]

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
    ]

    for variable in variables:
        try:
            stack_plot_qcd_cr(acc, outtag, region='sr_vbf_qcd_cr_detajj', variable=variable)
        except KeyError:
            print(f'Distribution not found: {variable}, skipping')

if __name__ == '__main__':
    main()
