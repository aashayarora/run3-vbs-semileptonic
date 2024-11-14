import ROOT as r
r.EnableImplicitMT()

import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use(hep.style.CMS)

import hist

class Plotter:
    def __init__(self, sig, bkg, data, xvar, xlabel, xbinning, yvar=None, ylabel=None, ybinning=None):
        self.sig = sig
        self.bkg = bkg
        self.data = data
        self.xvar = xvar
        self.xlabel = xlabel
        self.xbinning = xbinning
        
        self.yvar = yvar
        self.ylabel = ylabel
        self.ybinning = ybinning

    def plot1D(self):
        fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": (4, 1)})
        
        hep.histplot(self.bkg, ax=ax[0], histtype="fill", stack=True, label=bkg_samples_labels.values())
        hep.histplot(self.sig, label="Signal x 1000", ax=ax[0], histtype="step", color="red", linewidth=2, yerr=False)
        hep.histplot(self.data, label="Data", ax=ax[0], histtype="errorbar", color="black")
        hep.histplot(hist_ratio, color="black", ax=ax[1], histtype="errorbar")
        
        self.format_plot(True, ax[0])
        ax[0].legend()
        ax[0].set_xlabel("")
        ax[1].set_xlabel(self.xlabel)
        ax[0].set_ylabel("Events")
        ax[1].set_ylabel("Data / MC")
        ax[1].set_ylim(0.8, 1.2)
        ax[1].axhline(1, color="black", linestyle="--")

        plt.savefig(f"plots/{output_tag}/plot_{output_tag}_{histogram.output_file}.png")
        print("Saved", histogram.var)
        plt.close()
        


    def plot(self, output):
        if (self.yvar == None)
        
