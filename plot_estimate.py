import uproot
from matplotlib import pyplot as plt
import mplhep as hep

f = uproot.open("qcdestimate.root")
category = 'monojet'

for year in [2017, 2018]:
    nom = f['qcd_{category}_{year}']
    bin_up = f['qcd_monojet_2017_qcdbinning_monojet_2017Down']
    bin_dn = f['qcd_monojet_2017_qcdbinning_monojet_2017Up']
    fit_up = f['qcd_monojet_2017_qcdfit_monojet_2017Down']
    fit_dn = f['qcd_monojet_2017_qcdfit_monojet_2017Up']

    fig, ax = plt.subplots()
    hep.histplot(nom.values, nom.bins)
    fig.savefig(f"qcdestimate_{category}_{year}.pdf")