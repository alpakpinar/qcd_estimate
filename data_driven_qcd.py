#!/usr/bin/env python
import os
import pickle
import re
from collections import defaultdict

import numpy as np
import uproot
from bucoffea.plot.util import (URTH1, fig_ratio, klepto_load, merge_datasets,
                                merge_extensions, scale_xs_lumi)
from coffea import hist
from coffea.hist.export import export1d
from matplotlib import pyplot as plt

from fitlib import TFFit

colors = [
 'crimson',
 'darkorange',
 'navajowhite',
 'purple'
]

def exponential(x,a,b,c):
    ret = a * np.exp(-b*x) + c
    ret[ret<0] = 0
    return ret

def exponential2(x,a,b,c,d):
    ret = a * np.exp(-b*x -c*x**2) + d
    ret[ret<0] = 0
    return ret

from bucoffea.plot.style import matplotlib_rc

matplotlib_rc()

def ratio_unc(num, denom, dnum, ddenom):
    return np.hypot(
        dnum * (1/denom),
        ddenom * num / (denom*denom)
    )
data_err_opts = {
    # 'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'elinewidth': 1,
}

pjoin = os.path.join

def dd():
    return defaultdict(dict)


def make_templates(acc, outdir, tag, bins, region, dphi_cr=slice(0.0,0.5), dphi_sr=slice(0.5,None)):
    '''
    Creates the input templates for the fits.

    We start with a 2D distribution of mjj vs delta phi, which is just a coffea histogram.
    Then, we integrate dphi slices and different dataset combinations to get QCD MC, non-QCD MC and data
    in all the relevant regions.
    '''
    distribution = "mjj_vs_dphi_qcd"

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    

    parameters = defaultdict(dd)
    f = uproot.recreate(pjoin(outdir, f"templates_{region}_{tag}.root"))
    for year in [2017,2018]:
        # REBINNED
        h = acc[distribution].integrate("region", region)
        # TODO: Update the rebinning for mjj!
        h=h.rebin("mjj", hist.Bin('mjj',r'$M_{jj}$ (GeV)', bins))
        histos = {}
        histos["qcd"]    = h[re.compile(f"QCD.*HT.*{year}")].integrate("dataset")
        histos["nonqcd"] = h[re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')].integrate("dataset")
        histos["data"]   = h[re.compile(f"MET_{year}")].integrate("dataset")

        sr_sumw = {}
        sr_sumw2 = {}
        cr_sumw = {}
        cr_sumw2 = {}


        for name in ["qcd","nonqcd","data"]:
            sr_sumw[name],sr_sumw2[name] = histos[name].integrate("dphi",dphi_sr).values(sumw2=True,overflow="none")[()]
            cr_sumw[name],cr_sumw2[name] = histos[name].integrate("dphi",dphi_cr).values(sumw2=True,overflow="none")[()]
        
        tf = sr_sumw["qcd"]/cr_sumw["qcd"]
        dtf = ratio_unc(sr_sumw["qcd"], cr_sumw["qcd"], np.sqrt(sr_sumw2["qcd"]), np.sqrt(cr_sumw2["qcd"]))

        for name in ["qcd","nonqcd","data"]:
            f[f"{region}_{year}_cr_{name}"] = export1d(histos[name].integrate("dphi",dphi_cr))
            f[f"{region}_{year}_sr_{name}"] = export1d(histos[name].integrate("dphi",dphi_sr))
        f[f"{region}_{year}_tf"] = URTH1(
                    edges=histos["qcd"].axis("mjj").edges(),
                    sumw=np.r_[0,tf,0],
                    sumw2=np.r_[0,dtf**2,0]
                    )

def fit_tf(outdir, tag, region):
    '''
    Consume the input templates, make TFs and fit them.
    '''
    f = uproot.open(pjoin(outdir, f"templates_{region}_{tag}.root"))


    fits = {}
    for year in [2017,2018]:
        print("Fit",tag, year)
        h = f[f'{region}_{year}_tf']

        tf = h.values
        dtf = np.sqrt(h.variances)

        for i in range(len(dtf)):
            if dtf[i] == 0:
                dtf[i] = dtf[i-1]
                print(i, dtf[i])
        bins = h.allbins[1:-1]
        bins[-1,1] = bins[-1,0] + bins[-2,1] - bins[-2,0]
        dx = 0.5*np.diff(bins   , axis=1)[:,0]
        x  = 0.5*np.sum(bins, axis=1)
        
        guess = [0.5,1e-2,0]
        fits[year] = TFFit(
            x=x,
            y=tf,
            dy=dtf,
            fun=exponential,
            p0=guess

        )
        
        # Plot results
        fits[year].fit()
        fig, ax, rax = fig_ratio()
        ax.errorbar(
                    x,
                    tf,
                    xerr=dx,
                    yerr=dtf,
                    fmt="o",
                    label="Coarse-binned histogram",
                    color="k"
                    )
        
        xinterp = np.linspace(150, max(x), 100)
        nominal = fits[year].evaluate(xinterp, "best")
        ax.plot(
            xinterp,
            nominal,
            color='crimson',
            linestyle='-'
        )
        rax.errorbar(
            x,
            tf / fits[year].evaluate(x, "best"),
            dtf / fits[year].evaluate(x, "best"),
            fmt="o",
            color="k"
        )

        variations = fits[year].evaluate_all(xinterp)
        i = 0
        for var in set([re.sub("_(up|dn)","", x) for x in variations.keys()]):
            if 'best' in var:
                continue
            for direction in ['up','dn']:
                ax.fill_between(
                    xinterp,
                    nominal,
                    variations[f'{var}_{direction}'],
                    color=colors[i],
                    alpha=0.5,
                    label=var if direction=='up' else None
                )
                rax.fill_between(
                    xinterp,
                    nominal / nominal,
                    variations[f'{var}_{direction}'] / nominal,
                    color=colors[i],
                    alpha=0.5,
                    label=var if direction=='up' else None
                )
            i+=1
        # Aesthetics
        rax.set_ylim(0,4)
        rax.grid(linestyle="--")
        rax.set_ylabel("Histogram / fit")
        ax.set_ylabel("QCD MC transfer factor SR / CR")
        ax.set_xlabel(r"$M_{jj} \ (GeV)$")
        rax.set_xlabel(r"$M_{jj} \ (GeV)$")
        ax.legend()
        ax.set_yscale("log")
        ax.set_ylim(1e-5,1e0)
        ax.set_title(f"{tag}, {year}")
        fig.savefig(pjoin(outdir,f"tf_fit_{region}_{tag}_{year}.pdf"),bbox_inches='tight')
        plt.close(fig)

    with open(pjoin(outdir, f"tf_fit_{region}_{tag}.pkl"),"wb") as f:
        pickle.dump(fits, f)



def tf_variations(outdir, region):
    '''
    Nice plots of fit variations.
    '''
    x = np.linspace(200,5000,100)
    for year in [2017,2018]:
        fits = {}
        for file in  os.listdir(outdir):
            m = re.match(f"tf_fit_{region}_nominal_bin_([a-z,0-9]*).pkl",file)
            if not m:
                continue
            bintag = m.groups()[0]
            with open(pjoin(outdir, file),"rb") as f:
                fits[bintag] = pickle.load(f)[year]

        fig, ax, rax = fig_ratio()
        nominal = fits['nom'].evaluate(x,"best")
        for tag, fit in fits.items():
            ax.plot(
                     x,
                     fit.evaluate(x,"best"),
                     label=tag,
                     ls='-',
                     lw=2)
        
            rax.plot(
                x, 
                fit.evaluate(x,"best") / nominal,
                lw=2
            )

        env_dn, env_up = fits['nom'].envelope(x)
        rax.fill_between(
                    x,
                    env_dn / nominal,
                    env_up / nominal,
                    color="dodgerblue",
                    alpha=0.25,
                    label='Fit uncertainty')

        rax.set_ylim(0,2)
        rax.grid(linestyle="--")
        rax.set_ylabel("Variations / nominal")
        ax.set_ylabel("QCD MC transfer factor SR / CR")
        ax.set_xlabel(r"$M_{jj} \ (GeV)$")
        rax.set_xlabel(r"$M_{jj} \ (GeV)$")
        ax.legend()
        ax.set_yscale("log")
        ax.set_ylim(1e-4,1e1)
        ax.set_title(f"{year}")
        rax.legend()
        fig.savefig(pjoin(outdir,f"tf_variations_{region}_{year}.pdf"),bbox_inches='tight')
        plt.close(fig)

def histdiff(h1, h2):
    sumw = h1.values - h2.values
    sumw2 = h1.variances + h2.variances
    return sumw, sumw2
    

def tf_closure(outdir, region):
    '''
    Consumes the TFs and creates validation plots.
    '''
    x = np.linspace(250,1400,100)

    plotdir = pjoin(outdir, "closure")
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    for year in [2017,2018]:
        for cut in [0.2,0.3,0.4]:
            tag = f"closure_{cut}".replace('.','p')

            # Load fits
            fits = {}
            for file in  os.listdir(outdir):
                m = re.match(f"tf_fit_{region}_{tag}_bin_([a-z,0-9]*).pkl",file)
                if not m:
                    continue
                bintag = m.groups()[0]
                with open(pjoin(outdir, file),"rb") as f:
                    fits[bintag] = pickle.load(f)[year]
            # Load templates
            f = uproot.open(pjoin(outdir,f"templates_{region}_{tag}_bin_nom.root"))
            
            cr_qcd_sumw, cr_qcd_sumw2 = histdiff(f[f"{region}_{year}_cr_data"], f[f"{region}_{year}_cr_nonqcd"])
            sr_qcd_sumw, sr_qcd_sumw2 = histdiff(f[f"{region}_{year}_sr_data"], f[f"{region}_{year}_sr_nonqcd"])


            # sr_qcd_sumw  = sr_qcd_sumw - offset
            # print(year, cut, offset, sr_qcd_sumw)
            sr_qcd_mc_sumw = f[f"{region}_{year}_sr_qcd"].values
            sr_qcd_mc_sumw2 = f[f"{region}_{year}_sr_qcd"].variances

            bins = f[f"{region}_{year}_cr_data"].edges
            bins[-1,1] = bins[-1,0] + bins[-2,1] - bins[-2,0]
            dx = 0.5*np.diff(bins   , axis=1)
            x  = 0.5*np.sum(bins, axis=1)


            fig, ax, rax = fig_ratio()
            nominal = fits['nom'].evaluate(x,"best")
            ax.errorbar(
                x,
                y=sr_qcd_sumw[1:],
                yerr=np.sqrt(sr_qcd_sumw2[1:]),
                fmt='o',
                color='navy',
                label='Data - non-QCD in target'
            )
            ax.errorbar(
                x,
                y=cr_qcd_sumw[1:],
                yerr=np.sqrt(cr_qcd_sumw2[1:]),
                fmt='o',
                color='crimson',
                label='Data - non-QCD in reference'
            )
            ax.errorbar(
                x,
                y=sr_qcd_mc_sumw[1:],
                yerr=np.sqrt(sr_qcd_mc_sumw2[1:]),
                fmt='o',
                color='darkorange',
                label='QCD MC in target'
            )
            

            for fittag, fit in fits.items():
                ax.plot(
                        x,
                        fit.evaluate(x,"best") * cr_qcd_sumw[1:],
                        label=fittag,
                        ls='-',
                        lw=2)

            nominal = fits['nom'].evaluate(x,"best") * cr_qcd_sumw[1:]
            nominal_sumw2 = fits['nom'].evaluate(x,"best") * cr_qcd_sumw2[1:]
            env_dn, env_up = fits['nom'].envelope(x)  * cr_qcd_sumw[1:]

            ax.plot(
                x,
                nominal,
                label="Nominal prediction",
                ls='-',
                lw=2
            )
            ax.fill_between(
                x,
                env_dn,
                env_up,
                label="Prediction uncertainty",
                ls='-',
                lw=2
            )
            rax.errorbar(
                x,
                sr_qcd_sumw[1:] / nominal,
                yerr = np.sqrt(sr_qcd_sumw2[1:]) / nominal,
                fmt='o',
                color='navy'
            )
            rax.errorbar(
                x,
                sr_qcd_mc_sumw[1:] / nominal,
                yerr = np.sqrt(sr_qcd_mc_sumw2[1:]) / nominal,
                fmt='o',
                color='darkorange'
            )
            rax.fill_between(
                        x,
                        env_dn / nominal,
                        env_up / nominal,
                        color="dodgerblue",
                        alpha=0.25,
                        label='Fit uncertainty')

            rax.set_ylim(0,3)
            # rax.grid(linestyle="--")
            # rax.set_ylabel("Variations / nominal")
            # ax.set_ylabel("QCD MC transfer factor SR / CR")
            ax.set_xlabel(r"$M_{jj} \ (GeV)$")
            rax.set_xlabel(r"$M_{jj} \ (GeV)$")
            ax.legend()
            ax.set_yscale("log")
            ax.set_ylim(1e-4,1e8)
            # ax.set_title(f"{year}")
            # rax.legend()
            fig.savefig(pjoin(plotdir,f"tf_closure_{region}_{tag}_{year}.png"),bbox_inches='tight')
            plt.close(fig)

def tf_prediction(outdir,region):
    '''
    Consumes the fitted TFs and creates the final BG prediction.
    '''
    x = np.linspace(250,1400,100)

    plotdir = pjoin(outdir, "prediction")
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    fout = uproot.recreate(f"qcdestimate_{region}.root")
    for year in [2017,2018]:

        # Load fits
        fits = {}
        for file in  os.listdir(outdir):
            m = re.match(f"tf_fit_{region}_nominal_bin_([a-z,0-9]*).pkl",file)
            if not m:
                continue
            bintag = m.groups()[0]
            with open(pjoin(outdir, file),"rb") as f:
                fits[bintag] = pickle.load(f)[year]
        # Load templates
        f = uproot.open(pjoin(outdir,f"templates_{region}_nominal_bin_nom.root"))
        
        cr_qcd_sumw, cr_qcd_sumw2 = histdiff(f[f"{region}_{year}_cr_data"], f[f"{region}_{year}_cr_nonqcd"])

        sr_qcd_mc_sumw = f[f"{region}_{year}_sr_qcd"].values
        sr_qcd_mc_sumw2 = f[f"{region}_{year}_sr_qcd"].variances
        cr_qcd_mc_sumw = f[f"{region}_{year}_cr_qcd"].values
        cr_qcd_mc_sumw2 = f[f"{region}_{year}_cr_qcd"].variances

        bins = f[f"{region}_{year}_cr_data"].allbins[1:-1]
        bins[-1,1] = bins[-1,0] + bins[-2,1] - bins[-2,0]
        dx = 0.5*np.diff(bins   , axis=1)
        x  = 0.5*np.sum(bins, axis=1)

        fig, ax, rax = fig_ratio()
        nominal = fits['nom'].evaluate(x,"best") * cr_qcd_sumw
        nominal_sumw2 = fits['nom'].evaluate(x,"best") * cr_qcd_sumw2

        # Save nominal to file
        channel = 'vbf'
        fout[f'qcd_{channel}_{year}'] = (nominal, f[f"{region}_{year}_cr_data"].edges)

        # Plot QCD MC
        ax.errorbar(
            x,
            y=sr_qcd_mc_sumw,
            yerr=np.sqrt(sr_qcd_mc_sumw2),
            fmt='o',
            color='k',
            label='QCD MC in target'
        )
        
        # Plot binning variations
        for fittag, fit in fits.items():
            if fittag!='alt3':
                continue
            varied = fit.evaluate(x,"best") * cr_qcd_sumw
            # ax.fill_between(
                    # x,
                    # varied,
                    # 2*nominal -  varied,
                    # ls='-',
                    # color='dodgerblue',
                    # alpha=0.5,
                    # label='Binning uncertainty'
                    # )
            # Write binning variations to file
            fout[f'qcd_{channel}_{year}_qcdbinning_{channel}_{year}Up'] = URTH1(
                edges=np.unique(bins[:-1]),
                sumw=np.r_[0,varied],
                sumw2=np.r_[0,np.zeros(len(varied))],
            )
            fout[f'qcd_{channel}_{year}_qcdbinning_{channel}_{year}Down'] = URTH1(
                edges=np.unique(bins[:-1]),
                sumw=np.r_[0,(2*nominal-varied)],
                sumw2=np.r_[0,np.zeros(len(varied))],
            )

        ax.fill_between(
                    x,
                    nominal/1.25,
                    nominal*1.25,
                    ls='-',
                    color='crimson',
                    alpha=0.5,
                    label='Closure uncertainty'
                    )

        env_dn, env_up = fits['nom'].envelope(x)  * cr_qcd_sumw

        # Write fit variation envelopes to file
        fout[f'qcd_{channel}_{year}_qcdfit_{channel}_{year}Up'] = URTH1(
                edges=np.unique(bins[:-1]),
                sumw=np.r_[0,env_up],
                sumw2=np.r_[0,np.zeros(len(env_up))],
            )
        fout[f'qcd_{channel}_{year}_qcdfit_{channel}_{year}Down'] = URTH1(
                edges=np.unique(bins[:-1]),
                sumw=np.r_[0,env_dn],
                sumw2=np.r_[0,np.zeros(len(env_up))],
            )
        ax.plot(
            x,
            nominal,
            label="Nominal prediction",
            ls='-',
            lw=2
        )
        ax.fill_between(
            x,
            env_dn,
            env_up,
            label="Fit uncertainty",
            alpha=0.5,
            color="darkorange"
        )
        rax.errorbar(
            x,
            sr_qcd_mc_sumw / nominal,
            yerr = np.sqrt(sr_qcd_mc_sumw2) / nominal,
            fmt='o',
            color='k'
        )
        rax.fill_between(
                    x,
                    env_dn / nominal,
                    env_up / nominal,
                    color="darkorange",
                    alpha=0.25,
                    label='Fit uncertainty')

        for fittag, fit in fits.items():
            if fittag!='alt3':
                continue
            # rax.fill_between(
                    # x,
                    # fit.evaluate(x,"best") * cr_qcd_sumw / nominal,
                    # 2-fit.evaluate(x,"best") * cr_qcd_sumw / nominal,
                    # label=fittag,
                    # color='dodgerblue',
                    # alpha=0.25)

        rax.fill_between(
            x,
            1/1.25,
            1*1.25,
            ls='-',
            color='crimson',
            alpha=0.5,
            label='Closure uncertainty'
            )
        rax.set_ylim(0,3)
        ax.set_xlabel(r"$M_{jj} \ (GeV)$")
        rax.set_xlabel(r"$M_{jj} \ (GeV)$")
        rax.set_ylabel("Ratio to prediction")
        ax.legend()
        ax.set_yscale("log")
        ax.set_ylim(1e-4,1e8)
        fig.savefig(pjoin(plotdir,f"tf_prediction_{region}_{year}.pdf"),bbox_inches='tight')
        plt.close(fig)

def main():
    # Input handling
    indir = "./input/merged_2020-10-26_vbfhinv_03Sep20v7_qcd_estimation_very_loose_recoil_regions_detajj_cat"
    acc = klepto_load(indir)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')
    distributions = ["mjj_vs_dphi_qcd"]

    # Merging, scale, etc
    for distribution in distributions:
        acc.load(distribution)
        acc[distribution] = merge_extensions(acc[distribution], acc, reweight_pu=not ('nopu' in distribution))
        scale_xs_lumi(acc[distribution])
        acc[distribution] = merge_datasets(acc[distribution])
        acc[distribution].axis('dataset').sorting = 'integral'

    # Binnings for different regions, categorized by detajj
    bins = { 
        'sr_vbf_qcd_small_detajj' : {
             'nom' :  [200,400,600,900,1200,1500,2000,2750,3500,5000],
        },
        'sr_vbf_qcd_large_detajj' : {
             'nom' :  [200,400,600,900,1200,1500,2000,2750,3500,5000],
        },
    }

    outdir = pjoin('./output/',indir.split('/')[-1])
    
    # Estimate for each region is completely independent
    for region in ['sr_vbf_qcd_small_detajj', 'sr_vbf_qcd_large_detajj']:
        # Independent estimates also for for different bins
        for bintag, binvals in bins[region].items():
            tag =  f"nominal_bin_{bintag}"
            make_templates(
                            acc,
                            outdir,
                            tag, 
                            region=region, 
                            bins=binvals)
            fit_tf(
                    outdir,
                    tag, 
                    region)
            
            # For validation/closure testing, use variable delta phi cuts
            for cut in [0.2,0.3, 0.4]:
                tag = f"closure_{cut}_bin_{bintag}".replace('.','p')

                make_templates(
                                acc, 
                                outdir, 
                                tag, 
                                dphi_cr=slice(0.,cut), 
                                dphi_sr=slice(cut,0.5),
                                bins=binvals,
                                region=region
                                )
                fit_tf(
                        outdir, 
                        tag, 
                        region)

        tf_variations(outdir, region)
        # tf_closure(outdir, region)
        tf_prediction(outdir, region)
if __name__ == "__main__":
    main()
