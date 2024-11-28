import ROOT as r
import pandas as pd
import numpy as np
from scipy.stats import gamma
import os

from argparse import ArgumentParser

class DataCard:
    def __init__(self, data, sig, dnn, bdt, c2v, year=None, blind=True):
        self.blind = blind
        self.year = year
        self.c2v = c2v

        self.dnn = dnn
        self.bdt = bdt

        self.data_df = r.RDataFrame("Events", data).Filter("abs(weight) < 100")
        if (self.year != None):
            self.data_df = self.data_df.Filter(f"sample_year == \"{self.year}\"")
        self.data_pdf = pd.DataFrame(self.data_df.AsNumpy(["weight", "VBSBDTscore", "abcdnet_score"]))
        
        self.sig_base = sig
        self.sig_df = r.RDataFrame("Events", self.sig_base + "sig_MVA.root").Filter("abs(weight) < 100")
        if (self.year != None):
            self.sig_df = self.sig_df.Filter(f"sample_year == \"{self.year}\"")

        if (c2v != -1):
            self.sig_df = self.sig_df.Redefine("weight", f"weight * LHEReweightingWeight[{c2v}]")

        self.sig_pdf = pd.DataFrame(self.sig_df.AsNumpy(["weight", "VBSBDTscore", "abcdnet_score", "sample_year"]))
        
        self.total_sig = self.sig_pdf.weight.sum()

        self.sig_a = self.sig_pdf.query(f"VBSBDTscore > {bdt} & abcdnet_score > {dnn}")
        self.sig_b = self.sig_pdf.query(f"VBSBDTscore > {bdt} & abcdnet_score < {dnn}")
        self.sig_c = self.sig_pdf.query(f"VBSBDTscore < {bdt} & abcdnet_score > {dnn}")
        self.sig_d = self.sig_pdf.query(f"VBSBDTscore < {bdt} & abcdnet_score < {dnn}")

        self.uncertainties = {}
        self.datacard = ""

        self.data = self.get_data()
        self.sig = self.get_signal()

    def get_data(self):
        if self.blind:
            a = 1
            a_up = 1
            a_dn = 1
        else:
            a = round(self.data_pdf.query(f"VBSBDTscore > {self.bdt} & abcdnet_score > {self.dnn}").weight.sum())
            a_up = round(gamma.ppf(1- (1-0.9973) / 2, a+1))
            a_dn = round(gamma.ppf((1 - 0.9973) / 2, a)) if a > 0 else 0

        b = round(self.data_pdf.query(f"VBSBDTscore > {self.bdt} & abcdnet_score < {self.dnn}").weight.sum())
        b_up = round(gamma.ppf(1- (1-0.9973) / 2, b+1))
        b_dn = round(gamma.ppf((1 - 0.9973) / 2, b)) if b > 0 else 0
        
        c = round(self.data_pdf.query(f"VBSBDTscore < {self.bdt} & abcdnet_score > {self.dnn}").weight.sum())
        c_up = round(gamma.ppf(1- (1-0.9973) / 2, c+1))
        c_dn = round(gamma.ppf((1 - 0.9973) / 2, c)) if c > 0 else 0
        
        d = round(self.data_pdf.query(f"VBSBDTscore < {self.bdt} & abcdnet_score < {self.dnn}").weight.sum())
        d_up = round(gamma.ppf(1- (1-0.9973) / 2, d+1))
        d_dn = round(gamma.ppf((1 - 0.9973) / 2, d)) if d > 0 else 0

        return {"A": (a, a_dn, a_up), "B": (b, b_dn, b_up), "C": (c, c_dn, c_up), "D": (d, d_dn, d_up)}

    def get_signal(self):
        a = self.sig_a.weight.sum()
        b = self.sig_b.weight.sum()
        c = self.sig_c.weight.sum()
        d = self.sig_d.weight.sum()

        if a == 0 or b == 0 or c == 0 or d == 0:
            raise ValueError("Could not calculate signal yield: one of the regions has zero events")
        
        a_err = 1 + np.sqrt((self.sig_a.weight ** 2).sum()) / a
        b_err = 1 + np.sqrt((self.sig_b.weight ** 2).sum()) / b
        c_err = 1 + np.sqrt((self.sig_c.weight ** 2).sum()) / c
        d_err = 1 + np.sqrt((self.sig_d.weight ** 2).sum()) / d

        return {"A": (a, a_err), "B": (b, b_err), "C": (c, c_err), "D": (d, d_err)}


    def get_sf_variation(self, column):
        var_up = f"weight / {column}[0] * {column}[1]"
        var_dn = f"weight / {column}[0] * {column}[2]"
        pdf_var = pd.DataFrame(self.sig_df.Define("weight_up", var_up).Define("weight_dn", var_dn).AsNumpy(["weight_up", "weight_dn", "VBSBDTscore", "abcdnet_score"]))

        a = self.sig_a.weight.sum()
        b = self.sig_b.weight.sum()
        c = self.sig_c.weight.sum()
        d = self.sig_d.weight.sum()

        a_up = pdf_var.query(f"VBSBDTscore > {self.bdt} & abcdnet_score > {self.dnn}").weight_up.sum()
        b_up = pdf_var.query(f"VBSBDTscore > {self.bdt} & abcdnet_score < {self.dnn}").weight_up.sum()
        c_up = pdf_var.query(f"VBSBDTscore < {self.bdt} & abcdnet_score > {self.dnn}").weight_up.sum()
        d_up = pdf_var.query(f"VBSBDTscore < {self.bdt} & abcdnet_score < {self.dnn}").weight_up.sum()

        a_dn = pdf_var.query(f"VBSBDTscore > {self.bdt} & abcdnet_score > {self.dnn}").weight_dn.sum()
        b_dn = pdf_var.query(f"VBSBDTscore > {self.bdt} & abcdnet_score < {self.dnn}").weight_dn.sum()
        c_dn = pdf_var.query(f"VBSBDTscore < {self.bdt} & abcdnet_score > {self.dnn}").weight_dn.sum()
        d_dn = pdf_var.query(f"VBSBDTscore < {self.bdt} & abcdnet_score < {self.dnn}").weight_dn.sum()

        variations = []
        if a == 0 or b == 0 or c == 0 or d == 0:
            raise ValueError(f"Could not calculate uncertainty for {column}: one of the regions has zero events")
        
        variations_up = [abs(a - a_up) / a, abs(b - b_up) / b, abs(c - c_up) / c, abs(d - d_up) / d]
        variations_dn = [abs(a - a_dn) / a, abs(b - b_dn) / b, abs(c - c_dn) / c, abs(d - d_dn) / d]

        variations.append(variations_up)
        variations.append(variations_dn)
        
        def max_var(variations):
            deviations = [abs(1 - x) for x in variations]
            return variations[deviations.index(max(deviations))]
        
        return [(1 + abs(max_var(x))) for x in zip(*variations)]
    
    def add_sf_uncertainty(self, name, column):
        self.uncertainties[name] = self.get_sf_variation(column)

    def get_corr_variation(self, file):
        df_up = r.RDataFrame("Events", self.sig_base + "sig_" + file + "_up_MVA.root")
        df_dn = r.RDataFrame("Events", self.sig_base + "sig_" + file + "_down_MVA.root")
        
        if (self.c2v != -1):
            df_up = df_up.Redefine("weight", f"weight * LHEReweightingWeight[{self.c2v}]")
            df_dn = df_dn.Redefine("weight", f"weight * LHEReweightingWeight[{self.c2v}]")

        corr_df_up = pd.DataFrame(df_up.AsNumpy(["weight", "VBSBDTscore", "abcdnet_score"]))
        corr_df_dn = pd.DataFrame(df_dn.AsNumpy(["weight", "VBSBDTscore", "abcdnet_score"]))

        corr_df_up = corr_df_up.drop(corr_df_up.query("abs(weight) > 100").index)
        corr_df_dn = corr_df_dn.drop(corr_df_dn.query("abs(weight) > 100").index)

        a = self.sig_a.weight.sum()
        b = self.sig_b.weight.sum()
        c = self.sig_c.weight.sum()
        d = self.sig_d.weight.sum()

        a_up = corr_df_up.query(f"VBSBDTscore > {self.bdt} & abcdnet_score > {self.dnn}").weight.sum()
        b_up = corr_df_up.query(f"VBSBDTscore > {self.bdt} & abcdnet_score < {self.dnn}").weight.sum()
        c_up = corr_df_up.query(f"VBSBDTscore < {self.bdt} & abcdnet_score > {self.dnn}").weight.sum()
        d_up = corr_df_up.query(f"VBSBDTscore < {self.bdt} & abcdnet_score < {self.dnn}").weight.sum()

        a_dn = corr_df_dn.query(f"VBSBDTscore > {self.bdt} & abcdnet_score > {self.dnn}").weight.sum()
        b_dn = corr_df_dn.query(f"VBSBDTscore > {self.bdt} & abcdnet_score < {self.dnn}").weight.sum()
        c_dn = corr_df_dn.query(f"VBSBDTscore < {self.bdt} & abcdnet_score > {self.dnn}").weight.sum()
        d_dn = corr_df_dn.query(f"VBSBDTscore < {self.bdt} & abcdnet_score < {self.dnn}").weight.sum()

        if a == 0 or b == 0 or c == 0 or d == 0:
            raise ValueError(f"Could not calculate uncertainty for {file}: one of the regions has zero events")
        
        variations = []
        variations_up = [abs(a - a_up) / a, abs(b - b_up) / b, abs(c - c_up) / c, abs(d - d_up) / d]
        variations_dn = [abs(a - a_dn) / a, abs(b - b_dn) / b, abs(c - c_dn) / c, abs(d - d_dn) / d]

        variations.append(variations_up)
        variations.append(variations_dn)
        
        vars = [(1 + max(x)) for x in zip(*variations)]
        return vars

    def add_corr_uncertainty(self, name, file):
        self.uncertainties[name] = self.get_corr_variation(file)

    @staticmethod
    def format_header_line(name, *args):
        return ("{:<50}" + "{:<10}" + "{:<20}" * (len(args) - 1) + "\n").format(name, *args)

    @staticmethod
    def format_data_line(name, *args):
        return ("{:<50}" + "{:<10}" + "{:<20}" * (len(args) - 1) + "\n").format(name, *[f"{arg:.5f}" if isinstance(arg, float) else arg for arg in args])

    @staticmethod
    def format_uncertainty_line(name, *args):
        return ("{:<50}" + "{:<10}" + "{:<20}" * (len(args) - 1) + "\n").format(name, *[f"{arg:.5f}" if isinstance(arg, float) else arg for arg in args])

    @staticmethod
    def format_sig_uncertainty_line(name, *args):
        return ("{:<50}" + "{:<10}" + "{:<20}" * (len(args) + 4) + "\n").format(name, "lnN", "-", "-", "-", "-", *[f"{arg:.5f}" if isinstance(arg, float) else arg for arg in args])

    @staticmethod
    def format_rateparams_line(name, *args):
        return ("{:<25}" + "{:<20}" * (len(args) - 2) + "{}" + "\t" + "{}" + "\n").format(name, *args)

    def xbb_reweight(self, year):
        xbb_2016pre = ["CMS_vbsvvh1lep_bTagFitXbb_13TeV_16preVFP", f"{1 + (0.03010 * self.sig_a.query('sample_year == \"2016preVFP\"').weight.sum() / self.sig_a.weight.sum()):.5f}", f"{1 + (0.03010 * self.sig_b.query('sample_year == \"2016preVFP\"').weight.sum() / self.sig_b.weight.sum()):.5f}", f"{1 + (0.03010 * self.sig_c.query('sample_year == \"2016preVFP\"').weight.sum() / self.sig_c.weight.sum()):.5f}", f"{1 + (0.03010 * self.sig_d.query('sample_year == \"2016preVFP\"').weight.sum() / self.sig_d.weight.sum()):.5f}"]
        xbb_2016post = ["CMS_vbsvvh1lep_bTagFitXbb_13TeV_16postVFP", f"{1 + (0.01200 * self.sig_a.query('sample_year == \"2016postVFP\"').weight.sum() / self.sig_a.weight.sum()):.5f}", f"{1 + (0.01200 * self.sig_b.query('sample_year == \"2016postVFP\"').weight.sum() / self.sig_b.weight.sum()):.5f}", f"{1 + (0.01200 * self.sig_c.query('sample_year == \"2016postVFP\"').weight.sum() / self.sig_c.weight.sum()):.5f}", f"{1 + (0.01200 * self.sig_d.query('sample_year == \"2016postVFP\"').weight.sum() / self.sig_d.weight.sum()):.5f}"]
        xbb_2017 = ["CMS_vbsvvh1lep_bTagFitXbb_13TeV_17", f"{1 + (0.02560 * self.sig_a.query('sample_year == \"2017\"').weight.sum() / self.sig_a.weight.sum()):.5f}", f"{1 + (0.02560 * self.sig_b.query('sample_year == \"2017\"').weight.sum() / self.sig_b.weight.sum()):.5f}", f"{1 + (0.02560 * self.sig_c.query('sample_year == \"2017\"').weight.sum() / self.sig_c.weight.sum()):.5f}", f"{1 + (0.02560 * self.sig_d.query('sample_year == \"2017\"').weight.sum() / self.sig_d.weight.sum()):.5f}"]
        xbb_2018 = ["CMS_vbsvvh1lep_bTagFitXbb_13TeV_18", f"{1 + (0.07080 * self.sig_a.query('sample_year == \"2018\"').weight.sum() / self.sig_a.weight.sum()):.5f}", f"{1 + (0.07080 * self.sig_b.query('sample_year == \"2018\"').weight.sum() / self.sig_b.weight.sum()):.5f}", f"{1 + (0.07080 * self.sig_c.query('sample_year == \"2018\"').weight.sum() / self.sig_c.weight.sum()):.5f}", f"{1 + (0.07080 * self.sig_d.query('sample_year == \"2018\"').weight.sum() / self.sig_d.weight.sum()):.5f}"]
        if (year == "2016preVFP"):
            return (xbb_2016pre)
        elif (year == "2016postVFP"):
            return (xbb_2016post)
        elif (year == "2017"):
            return (xbb_2017)
        elif (year == "2018"):
            return (xbb_2018)
        else:
            return (xbb_2016pre, xbb_2016post, xbb_2017, xbb_2018)

    def generate_datacard(self):
        self.datacard += "imax 4 number of channels\n"
        self.datacard += "jmax 1 number of backgrounds\n"
        self.datacard += "kmax * number of nuisance parameters\n"
        self.datacard += "-" * 150 + "\n"
        self.datacard += self.format_header_line("bin", "", "A", "B", "C", "D")
        self.datacard += self.format_data_line("observation", "", self.data["A"][0], self.data["B"][0], self.data["C"][0], self.data["D"][0])
        self.datacard += "-" * 150 + "\n"
        self.datacard += self.format_header_line("bin", "", "A", "B", "C", "D", "A", "B", "C", "D")
        self.datacard += self.format_header_line("process", "", "TotalBkg_OneLep", "TotalBkg_OneLep", "TotalBkg_OneLep", "TotalBkg_OneLep", "TotalSig", "TotalSig", "TotalSig", "TotalSig")
        self.datacard += self.format_header_line("process", "", "1", "1", "1", "1", "0", "0", "0", "0")
        self.datacard += self.format_data_line("rate", "", "1", "1", "1", "1", self.sig["A"][0], self.sig["B"][0], self.sig["C"][0], self.sig["D"][0])
        self.datacard += "-" * 150 + "\n"
        self.datacard += self.format_uncertainty_line("CMS_vbsvvh1lep_control_abcd_syst", "lnN", "1.35", "-", "-", "-", "-", "-", "-", "-")
        self.datacard += self.format_uncertainty_line("lumi_13TeV_correlated", "lnN", "-", "-", "-", "-", "1.016", "1.016", "1.016", "1.016")
        self.datacard += self.format_uncertainty_line("CMS_vbsvvh1lep_signal_RegionA", "lnN", "-", "-", "-", "-", self.sig["A"][1], "-", "-", "-")
        self.datacard += self.format_uncertainty_line("CMS_vbsvvh1lep_signal_RegionB", "lnN", "-", "-", "-", "-", "-", self.sig["B"][1], "-", "-")
        self.datacard += self.format_uncertainty_line("CMS_vbsvvh1lep_signal_RegionC", "lnN", "-", "-", "-", "-", "-", "-", self.sig["C"][1], "-")
        self.datacard += self.format_uncertainty_line("CMS_vbsvvh1lep_signal_RegionD", "lnN", "-", "-", "-", "-", "-", "-", "-", self.sig["D"][1])
        for name, values in self.uncertainties.items():
            self.datacard += self.format_sig_uncertainty_line(name, *values)
        xbb_weights = self.xbb_reweight(self.year)
        for xbb in xbb_weights:
            self.datacard += self.format_sig_uncertainty_line(*xbb)
        self.datacard += "-" * 150 + "\n"
        self.datacard += self.format_rateparams_line("A_OneLep rateParam", "A", "TotalBkg_OneLep", "(@0*@1/@2)", "B_OneLep,C_OneLep,D_OneLep")
        self.datacard += self.format_rateparams_line("B_OneLep rateParam", "B", "TotalBkg_OneLep", self.data["B"][0], f"[{self.data["B"][1]},{self.data["B"][2]}]")
        self.datacard += self.format_rateparams_line("C_OneLep rateParam", "C", "TotalBkg_OneLep", self.data["C"][0], f"[{self.data["C"][1]},{self.data["C"][2]}]")
        self.datacard += self.format_rateparams_line("D_OneLep rateParam", "D", "TotalBkg_OneLep", self.data["D"][0], f"[{self.data["D"][1]},{self.data["D"][2]}]")

    def __str__(self):
        return self.datacard

if __name__ == "__main__":
    argparser = ArgumentParser()
    argparser.add_argument("--data", type=str, required=True, help="Path to data file")
    argparser.add_argument("--sig", type=str, required=True, help="Path to signal file")
    argparser.add_argument("--dnn", type=float, required=True, help="DNN cut")
    argparser.add_argument("--bdt", type=float, required=True, help="BDT cut")
    argparser.add_argument("--c2v", type=int, required=False, help="c2v value")
    argparser.add_argument("--year", type=int, required=False, help="year")
    argparser.add_argument("--output_dir", type=str, required=False, help="output directory")
    argparser.add_argument("--wzscan", action="store_true", help="add WZ scan uncertainties")
    args = argparser.parse_args()
    
    datacard = DataCard(data=args.data, sig=args.sig, dnn=args.dnn, bdt=args.bdt, c2v=args.c2v, year=args.year)
    datacard.add_sf_uncertainty("CMS_LHE_weights_pdf_vbsvvh", "LHEWeights_pdf")
    datacard.add_sf_uncertainty("CMS_vbsvvh_puWeight", "pileup_weight")
    datacard.add_sf_uncertainty("CMS_ttH_elec_reco", "electron_scale_factors_Reco")
    datacard.add_sf_uncertainty("CMS_ttH_elec_recotoloose", "electron_scale_factors_ID")
    datacard.add_sf_uncertainty("CMS_ttH_elec_trig", "electron_scale_factors_trigger")
    datacard.add_sf_uncertainty("CMS_ttH_elec_loosetoiso", "electron_scale_factors_ttHISO")
    datacard.add_sf_uncertainty("CMS_ttH_elec_isototight", "electron_scale_factors_ttHID")
    datacard.add_sf_uncertainty("CMS_ttH_muon_recotoloose", "muon_scale_factors_ID")
    datacard.add_sf_uncertainty("CMS_ttH_muon_trig", "muon_scale_factors_trigger")
    datacard.add_sf_uncertainty("CMS_ttH_muon_loosetoiso", "muon_scale_factors_ttHISO")
    datacard.add_sf_uncertainty("CMS_ttH_muon_isototight", "muon_scale_factors_ttHID")
    datacard.add_sf_uncertainty("CMS_PSWeight_FSR_vbsvvh", "PSWeight_FSR")
    datacard.add_sf_uncertainty("CMS_PSWeight_ISR_vbsvvh", "PSWeight_ISR")
    datacard.add_sf_uncertainty("CMS_PrefireWeight_13TeV", "L1PreFiringWeight")
    datacard.add_sf_uncertainty("CMS_vbsvvh_puJetID", "pileupid_weight")
    datacard.add_corr_uncertainty("CMS_scale_j_Absolute_13TeV", "jec_absolute")
    datacard.add_corr_uncertainty("CMS_scale_j_Absolute_2016postVFP_13TeV", "jec_absolute_2016postVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_Absolute_2016preVFP_13TeV", "jec_absolute_2016preVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_Absolute_2017_13TeV", "jec_absolute_2017")
    datacard.add_corr_uncertainty("CMS_scale_j_Absolute_2018_13TeV", "jec_absolute_2018")
    datacard.add_corr_uncertainty("CMS_scale_j_BBEC1_13TeV", "jec_bbec1")
    datacard.add_corr_uncertainty("CMS_scale_j_BBEC1_2016postVFP_13TeV", "jec_bbec1_2016postVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_BBEC1_2016preVFP_13TeV", "jec_bbec1_2016preVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_BBEC1_2017_13TeV", "jec_bbec1_2017")
    datacard.add_corr_uncertainty("CMS_scale_j_BBEC1_2018_13TeV", "jec_bbec1_2018")
    datacard.add_corr_uncertainty("CMS_scale_j_EC2_13TeV", "jec_ec2")
    datacard.add_corr_uncertainty("CMS_scale_j_EC2_2016postVFP_13TeV", "jec_ec2_2016postVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_EC2_2016preVFP_13TeV", "jec_ec2_2016preVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_EC2_2017_13TeV", "jec_ec2_2017")
    datacard.add_corr_uncertainty("CMS_scale_j_EC2_2018_13TeV", "jec_ec2_2018")
    datacard.add_corr_uncertainty("CMS_scale_j_FlavorQCD_13TeV", "jec_flavorqcd")
    datacard.add_corr_uncertainty("CMS_scale_j_HF_13TeV", "jec_hf")
    datacard.add_corr_uncertainty("CMS_scale_j_HF_2016postVFP_13TeV", "jec_hf_2016postVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_HF_2016preVFP_13TeV", "jec_hf_2016preVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_HF_2017_13TeV", "jec_hf_2017")
    datacard.add_corr_uncertainty("CMS_scale_j_HF_2018_13TeV", "jec_hf_2018")
    datacard.add_corr_uncertainty("CMS_scale_j_RelativeBal_13TeV", "jec_relativebal")
    datacard.add_corr_uncertainty("CMS_scale_j_RelativeSample_2016postVFP_13TeV", "jec_relativesample_2016postVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_RelativeSample_2016preVFP_13TeV", "jec_relativesample_2016preVFP")
    datacard.add_corr_uncertainty("CMS_scale_j_RelativeSample_2017_13TeV", "jec_relativesample_2017")
    datacard.add_corr_uncertainty("CMS_scale_j_RelativeSample_2018_13TeV", "jec_relativesample_2018")
    datacard.add_corr_uncertainty("CMS_res_j_13TeV", "jer")
    datacard.add_corr_uncertainty("CMS_metUncl_13TeV", "met_unclustered")
    datacard.add_corr_uncertainty("CMS_jms_pnetreg", "jms")
    datacard.add_corr_uncertainty("CMS_jmr_pnetreg", "jmr")
    datacard.add_sf_uncertainty("CMS_LHE_weights_scale_muF_vbsvvh", "LHEScaleWeight_muF")
    datacard.add_sf_uncertainty("CMS_LHE_weights_scale_muR_vbsvvh", "LHEScaleWeight_muR")
    datacard.add_sf_uncertainty("CMS_btagWeightDeepJet_HF_13Tev", "btagging_scale_factors_HF")
    datacard.add_sf_uncertainty("CMS_btagWeightDeepJet_LF_13Tev", "btagging_scale_factors_LF")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_qTagWeightXWqq_13TeV_16preVFP", "particlenet_w_weight_2016preVFP")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_qTagWeightXWqq_13TeV_16postVFP", "particlenet_w_weight_2016postVFP")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_qTagWeightXWqq_13TeV_17", "particlenet_w_weight_2017")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_qTagWeightXWqq_13TeV_18", "particlenet_w_weight_2018")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_bTagWeightXbb_13TeV_16preVFP", "particlenet_h_weight_2016preVFP")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_bTagWeightXbb_13TeV_16postVFP", "particlenet_h_weight_2016postVFP")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_bTagWeightXbb_13TeV_17", "particlenet_h_weight_2017")
    datacard.add_sf_uncertainty("CMS_vbsvvh1lep_bTagWeightXbb_13TeV_18", "particlenet_h_weight_2018")
    datacard.generate_datacard()

    if args.wzscan:
        names = [ 
            "scan_C2W_m2p0_C2Z_m2p0", 
            "scan_C2W_m2p0_C2Z_m1p0", 
            "scan_C2W_m2p0_C2Z_m0p5", 
            "scan_C2W_m2p0_C2Z_0p0", 
            "scan_C2W_m2p0_C2Z_0p5", 
            "scan_C2W_m2p0_C2Z_1p0", 
            "scan_C2W_m2p0_C2Z_1p5", 
            "scan_C2W_m2p0_C2Z_2p0", 
            "scan_C2W_m2p0_C2Z_2p5", 
            "scan_C2W_m2p0_C2Z_3p0", 
            "scan_C2W_m2p0_C2Z_4p0", 
            "scan_C2W_m1p0_C2Z_m2p0", 
            "scan_C2W_m1p0_C2Z_m1p0", 
            "scan_C2W_m1p0_C2Z_m0p5", 
            "scan_C2W_m1p0_C2Z_0p0", 
            "scan_C2W_m1p0_C2Z_0p5", 
            "scan_C2W_m1p0_C2Z_1p0", 
            "scan_C2W_m1p0_C2Z_1p5", 
            "scan_C2W_m1p0_C2Z_2p0", 
            "scan_C2W_m1p0_C2Z_2p5", 
            "scan_C2W_m1p0_C2Z_3p0", 
            "scan_C2W_m1p0_C2Z_4p0", 
            "scan_C2W_m0p5_C2Z_m2p0", 
            "scan_C2W_m0p5_C2Z_m1p0", 
            "scan_C2W_m0p5_C2Z_m0p5", 
            "scan_C2W_m0p5_C2Z_0p0", 
            "scan_C2W_m0p5_C2Z_0p5", 
            "scan_C2W_m0p5_C2Z_1p0", 
            "scan_C2W_m0p5_C2Z_1p5", 
            "scan_C2W_m0p5_C2Z_2p0", 
            "scan_C2W_m0p5_C2Z_2p5", 
            "scan_C2W_m0p5_C2Z_3p0", 
            "scan_C2W_m0p5_C2Z_4p0", 
            "scan_C2W_0p0_C2Z_m2p0", 
            "scan_C2W_0p0_C2Z_m1p0", 
            "scan_C2W_0p0_C2Z_m0p5", 
            "scan_C2W_0p0_C2Z_0p0", 
            "scan_C2W_0p0_C2Z_0p5", 
            "scan_C2W_0p0_C2Z_1p0", 
            "scan_C2W_0p0_C2Z_1p5", 
            "scan_C2W_0p0_C2Z_2p0", 
            "scan_C2W_0p0_C2Z_2p5", 
            "scan_C2W_0p0_C2Z_3p0", 
            "scan_C2W_0p0_C2Z_4p0", 
            "scan_C2W_0p5_C2Z_m2p0", 
            "scan_C2W_0p5_C2Z_m1p0", 
            "scan_C2W_0p5_C2Z_m0p5", 
            "scan_C2W_0p5_C2Z_0p0", 
            "scan_C2W_0p5_C2Z_0p5", 
            "scan_C2W_0p5_C2Z_1p0", 
            "scan_C2W_0p5_C2Z_1p5", 
            "scan_C2W_0p5_C2Z_2p0", 
            "scan_C2W_0p5_C2Z_2p5", 
            "scan_C2W_0p5_C2Z_3p0", 
            "scan_C2W_0p5_C2Z_4p0", 
            "scan_C2W_1p0_C2Z_m2p0", 
            "scan_C2W_1p0_C2Z_m1p0", 
            "scan_C2W_1p0_C2Z_m0p5", 
            "scan_C2W_1p0_C2Z_0p0", 
            "scan_C2W_1p0_C2Z_0p5", 
            "scan_C2W_1p0_C2Z_1p0", 
            "scan_C2W_1p0_C2Z_1p5", 
            "scan_C2W_1p0_C2Z_2p0", 
            "scan_C2W_1p0_C2Z_2p5", 
            "scan_C2W_1p0_C2Z_3p0", 
            "scan_C2W_1p0_C2Z_4p0", 
            "scan_C2W_1p5_C2Z_m2p0", 
            "scan_C2W_1p5_C2Z_m1p0", 
            "scan_C2W_1p5_C2Z_m0p5", 
            "scan_C2W_1p5_C2Z_0p0", 
            "scan_C2W_1p5_C2Z_0p5", 
            "scan_C2W_1p5_C2Z_1p0", 
            "scan_C2W_1p5_C2Z_1p5", 
            "scan_C2W_1p5_C2Z_2p0", 
            "scan_C2W_1p5_C2Z_2p5", 
            "scan_C2W_1p5_C2Z_3p0", 
            "scan_C2W_1p5_C2Z_4p0", 
            "scan_C2W_2p0_C2Z_m2p0", 
            "scan_C2W_2p0_C2Z_m1p0", 
            "scan_C2W_2p0_C2Z_m0p5", 
            "scan_C2W_2p0_C2Z_0p0", 
            "scan_C2W_2p0_C2Z_0p5", 
            "scan_C2W_2p0_C2Z_1p0", 
            "scan_C2W_2p0_C2Z_1p5", 
            "scan_C2W_2p0_C2Z_2p0", 
            "scan_C2W_2p0_C2Z_2p5", 
            "scan_C2W_2p0_C2Z_3p0", 
            "scan_C2W_2p0_C2Z_4p0", 
            "scan_C2W_2p5_C2Z_m2p0", 
            "scan_C2W_2p5_C2Z_m1p0", 
            "scan_C2W_2p5_C2Z_m0p5", 
            "scan_C2W_2p5_C2Z_0p0", 
            "scan_C2W_2p5_C2Z_0p5", 
            "scan_C2W_2p5_C2Z_1p0", 
            "scan_C2W_2p5_C2Z_1p5", 
            "scan_C2W_2p5_C2Z_2p0", 
            "scan_C2W_2p5_C2Z_2p5", 
            "scan_C2W_2p5_C2Z_3p0", 
            "scan_C2W_2p5_C2Z_4p0", 
            "scan_C2W_3p0_C2Z_m2p0", 
            "scan_C2W_3p0_C2Z_m1p0", 
            "scan_C2W_3p0_C2Z_m0p5", 
            "scan_C2W_3p0_C2Z_0p0", 
            "scan_C2W_3p0_C2Z_0p5", 
            "scan_C2W_3p0_C2Z_1p0", 
            "scan_C2W_3p0_C2Z_1p5", 
            "scan_C2W_3p0_C2Z_2p0", 
            "scan_C2W_3p0_C2Z_2p5", 
            "scan_C2W_3p0_C2Z_3p0", 
            "scan_C2W_3p0_C2Z_4p0", 
            "scan_C2W_4p0_C2Z_m2p0", 
            "scan_C2W_4p0_C2Z_m1p0", 
            "scan_C2W_4p0_C2Z_m0p5", 
            "scan_C2W_4p0_C2Z_0p0", 
            "scan_C2W_4p0_C2Z_0p5", 
            "scan_C2W_4p0_C2Z_1p0", 
            "scan_C2W_4p0_C2Z_1p5", 
            "scan_C2W_4p0_C2Z_2p0", 
            "scan_C2W_4p0_C2Z_2p5", 
            "scan_C2W_4p0_C2Z_3p0", 
            "scan_C2W_4p0_C2Z_4p0",
        ]

        myrange=list(range(0, len(names)-1))
        # myrange.insert(names.index("scan_CV_1p0_C2V_2p0_C3_1p0"), -1)
        myrange.insert(names.index("scan_C2W_2p0_C2Z_2p0"), -1)
        datacard_name = names[myrange.index(args.c2v)]

    else:
        names = [
            "scan_CV_1p0_C2V_m2p0_C3_1p0",
            "scan_CV_1p0_C2V_m1p75_C3_1p0",
            "scan_CV_1p0_C2V_m1p5_C3_1p0",
            "scan_CV_1p0_C2V_m1p25_C3_1p0",
            "scan_CV_1p0_C2V_m1p0_C3_1p0",
            "scan_CV_1p0_C2V_m0p75_C3_1p0",
            "scan_CV_1p0_C2V_m0p5_C3_1p0",
            "scan_CV_1p0_C2V_m0p25_C3_1p0",
            "scan_CV_1p0_C2V_0p0_C3_1p0",
            "scan_CV_1p0_C2V_0p1_C3_1p0",
            "scan_CV_1p0_C2V_0p2_C3_1p0",
            "scan_CV_1p0_C2V_0p3_C3_1p0",
            "scan_CV_1p0_C2V_0p4_C3_1p0",
            "scan_CV_1p0_C2V_0p5_C3_1p0",
            "scan_CV_1p0_C2V_0p6_C3_1p0",
            "scan_CV_1p0_C2V_0p7_C3_1p0",
            "scan_CV_1p0_C2V_0p8_C3_1p0",
            "scan_CV_1p0_C2V_0p9_C3_1p0",
            "scan_CV_1p0_C2V_1p0_C3_1p0",
            "scan_CV_1p0_C2V_1p1_C3_1p0",
            "scan_CV_1p0_C2V_1p2_C3_1p0",
            "scan_CV_1p0_C2V_1p3_C3_1p0",
            "scan_CV_1p0_C2V_1p4_C3_1p0",
            "scan_CV_1p0_C2V_1p5_C3_1p0",
            "scan_CV_1p0_C2V_1p6_C3_1p0",
            "scan_CV_1p0_C2V_1p7_C3_1p0",
            "scan_CV_1p0_C2V_1p8_C3_1p0",
            "scan_CV_1p0_C2V_1p9_C3_1p0",
            "scan_CV_1p0_C2V_2p0_C3_1p0",
            "scan_CV_1p0_C2V_2p25_C3_1p0",
            "scan_CV_1p0_C2V_2p5_C3_1p0",
            "scan_CV_1p0_C2V_2p75_C3_1p0",
            "scan_CV_1p0_C2V_3p0_C3_1p0",
            "scan_CV_1p0_C2V_3p25_C3_1p0",
            "scan_CV_1p0_C2V_3p5_C3_1p0",
            "scan_CV_1p0_C2V_3p75_C3_1p0",
            "scan_CV_1p0_C2V_4p0_C3_1p0"
        ]

        myrange=list(range(0, len(names)-1))
        myrange.insert(names.index("scan_CV_1p0_C2V_2p0_C3_1p0"), -1)
        datacard_name = names[myrange.index(args.c2v)]

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(f"{args.output_dir}/{datacard_name}.dat" , "w") as f:
        f.write(str(datacard))