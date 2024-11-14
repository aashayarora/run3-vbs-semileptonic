#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"

#include <vector>
#include <iostream>

using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;
using ROOT::RDF::RSampleInfo;

float looseDFBtagWP(std::string year){
    if(year == "2016preVFP")
        return 0.0508;
    if(year == "2016postVFP")
        return 0.0480;
    if(year == "2017")
        return 0.0532;
    if(year == "2018")
        return 0.0490;
    return -1;
}

float mediumDFBtagWP(std::string year){
    if(year == "2016preVFP")
        return 0.2598;
    if(year == "2016postVFP")
        return 0.2489;
    if(year == "2017")
        return 0.3040;
    if(year == "2018")
        return 0.2783;
    return -1;
}

float tightDFBtagWP(std::string year){
    if(year == "2016preVFP")
        return 0.6502;
    if(year == "2016postVFP")
        return 0.6377;
    if(year == "2017")
        return 0.7476;
    if(year == "2018")
        return 0.7100;
    return -1;
}

int main(int arc, char** argv){
    ROOT::EnableImplicitMT(64);
    ROOT::RDataFrame df = ROOT::RDF::Experimental::FromSpec("../input/sig.json");
    ROOT::RDF::Experimental::AddProgressBar(df);

    auto df0 = df.DefinePerSample("xsec_weight", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("xsec_weight");})
            .DefinePerSample("sample_category", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_category");})
            .DefinePerSample("sample_type", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_type");})
            .DefinePerSample("sample_year", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_year");})
            .Define("isData", "sample_category == \"data\"")
            .Define("is2016", "sample_year == \"2016preVFP\" || sample_year == \"2016postVFP\"")
            .Define("is2017", "sample_year == \"2017\"")
            .Define("is2018", "sample_year == \"2018\"");

    auto df1 = df0.Define("ak4tightBjetScore", tightDFBtagWP, {"sample_year"})
                .Define("ak4mediumBjetScore", mediumDFBtagWP, {"sample_year"})
                .Define("ak4looseBjetScore", looseDFBtagWP, {"sample_year"})
                .Define("tightBJets", "Jet_btagDeepFlavB > ak4tightBjetScore && Jet_hadronFlavour == 5") 
                .Define("tightbjet_pt", "Jet_pt[tightBJets]")
                .Define("tightbjet_eta", "Jet_eta[tightBJets]")
                .Define("looseBJets", "Jet_btagDeepFlavB > ak4looseBjetScore && Jet_btagDeepFlavB < ak4tightBjetScore && Jet_hadronFlavour == 5")
                .Define("loosebjet_pt", "Jet_pt[looseBJets]")
                .Define("loosebjet_eta", "Jet_eta[looseBJets]")
                .Define("notagBJets", "Jet_btagDeepFlavB < ak4looseBjetScore && Jet_hadronFlavour == 5")
                .Define("notagbjet_pt", "Jet_pt[notagBJets]")
                .Define("notagbjet_eta", "Jet_eta[notagBJets]")
                .Define("tightCJets", "Jet_btagDeepFlavB > ak4tightBjetScore && Jet_hadronFlavour == 4")
                .Define("tightcjet_pt", "Jet_pt[tightCJets]")
                .Define("tightcjet_eta", "Jet_eta[tightCJets]")
                .Define("looseCJets", "Jet_btagDeepFlavB > ak4looseBjetScore && Jet_btagDeepFlavB < ak4tightBjetScore && Jet_hadronFlavour == 4")
                .Define("loosecjet_pt", "Jet_pt[looseCJets]")
                .Define("loosecjet_eta", "Jet_eta[looseCJets]")
                .Define("notagCJets", "Jet_btagDeepFlavB < ak4looseBjetScore && Jet_hadronFlavour == 4")
                .Define("notagcjet_pt", "Jet_pt[notagCJets]")
                .Define("notagcjet_eta", "Jet_eta[notagCJets]")
                .Define("tightLJets", "Jet_btagDeepFlavB > ak4tightBjetScore && Jet_hadronFlavour == 0")
                .Define("tightljet_pt", "Jet_pt[tightLJets]")
                .Define("tightljet_eta", "Jet_eta[tightLJets]")
                .Define("looseLJets", "Jet_btagDeepFlavB > ak4looseBjetScore && Jet_btagDeepFlavB < ak4tightBjetScore && Jet_hadronFlavour == 0")
                .Define("looseljet_pt", "Jet_pt[looseLJets]")
                .Define("looseljet_eta", "Jet_eta[looseLJets]")
                .Define("notagLJets", "Jet_btagDeepFlavB < ak4looseBjetScore && Jet_hadronFlavour == 0")
                .Define("notagljet_pt", "Jet_pt[notagLJets]")
                .Define("notagljet_eta", "Jet_eta[notagLJets]");

    double pt_bins[11] = {15, 20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
    double eta_bins[5] = {-2.5, -0.9, 0, 0.9, 2.5};

    auto hist_btight_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_btight_2016preVFP", "hist_btight_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightbjet_pt", "tightbjet_eta");    
    auto hist_bloose_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_bloose_2016preVFP", "hist_bloose_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosebjet_pt", "loosebjet_eta");
    auto hist_bnotag_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_bnotag_2016preVFP", "hist_bnotag_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagbjet_pt", "notagbjet_eta");

    auto hist_ctight_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_ctight_2016preVFP", "hist_ctight_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightcjet_pt", "tightcjet_eta");
    auto hist_cloose_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_cloose_2016preVFP", "hist_cloose_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosecjet_pt", "loosecjet_eta");
    auto hist_cnotag_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_cnotag_2016preVFP", "hist_cnotag_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagcjet_pt", "notagcjet_eta");

    auto hist_ltight_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_ltight_2016preVFP", "hist_ltight_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightljet_pt", "tightljet_eta");
    auto hist_lloose_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_lloose_2016preVFP", "hist_lloose_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "looseljet_pt", "looseljet_eta");
    auto hist_lnotag_2016preVFP = df1.Filter("sample_year == \"2016preVFP\"").Histo2D({"hist_lnotag_2016preVFP", "hist_lnotag_2016preVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagljet_pt", "notagljet_eta");

    auto hist_btight_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_btight_2016postVFP", "hist_btight_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightbjet_pt", "tightbjet_eta");    
    auto hist_bloose_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_bloose_2016postVFP", "hist_bloose_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosebjet_pt", "loosebjet_eta");
    auto hist_bnotag_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_bnotag_2016postVFP", "hist_bnotag_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagbjet_pt", "notagbjet_eta");

    auto hist_ctight_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_ctight_2016postVFP", "hist_ctight_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightcjet_pt", "tightcjet_eta");
    auto hist_cloose_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_cloose_2016postVFP", "hist_cloose_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosecjet_pt", "loosecjet_eta");
    auto hist_cnotag_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_cnotag_2016postVFP", "hist_cnotag_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagcjet_pt", "notagcjet_eta");

    auto hist_ltight_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_ltight_2016postVFP", "hist_ltight_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightljet_pt", "tightljet_eta");
    auto hist_lloose_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_lloose_2016postVFP", "hist_lloose_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "looseljet_pt", "looseljet_eta");
    auto hist_lnotag_2016postVFP = df1.Filter("sample_year == \"2016postVFP\"").Histo2D({"hist_lnotag_2016postVFP", "hist_lnotag_2016postVFP;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagljet_pt", "notagljet_eta");
    
    auto hist_btight_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_btight_2017", "hist_btight_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightbjet_pt", "tightbjet_eta");    
    auto hist_bloose_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_bloose_2017", "hist_bloose_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosebjet_pt", "loosebjet_eta");
    auto hist_bnotag_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_bnotag_2017", "hist_bnotag_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagbjet_pt", "notagbjet_eta");

    auto hist_ctight_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_ctight_2017", "hist_ctight_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightcjet_pt", "tightcjet_eta");
    auto hist_cloose_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_cloose_2017", "hist_cloose_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosecjet_pt", "loosecjet_eta");
    auto hist_cnotag_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_cnotag_2017", "hist_cnotag_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagcjet_pt", "notagcjet_eta");

    auto hist_ltight_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_ltight_2017", "hist_ltight_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightljet_pt", "tightljet_eta");
    auto hist_lloose_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_lloose_2017", "hist_lloose_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "looseljet_pt", "looseljet_eta");
    auto hist_lnotag_2017 = df1.Filter("sample_year == \"2017\"").Histo2D({"hist_lnotag_2017", "hist_lnotag_2017;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagljet_pt", "notagljet_eta");
   
    auto hist_btight_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_btight_2018", "hist_btight_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightbjet_pt", "tightbjet_eta");    
    auto hist_bloose_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_bloose_2018", "hist_bloose_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosebjet_pt", "loosebjet_eta");
    auto hist_bnotag_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_bnotag_2018", "hist_bnotag_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagbjet_pt", "notagbjet_eta");

    auto hist_ctight_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_ctight_2018", "hist_ctight_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightcjet_pt", "tightcjet_eta");
    auto hist_cloose_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_cloose_2018", "hist_cloose_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "loosecjet_pt", "loosecjet_eta");
    auto hist_cnotag_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_cnotag_2018", "hist_cnotag_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagcjet_pt", "notagcjet_eta");

    auto hist_ltight_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_ltight_2018", "hist_ltight_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "tightljet_pt", "tightljet_eta");
    auto hist_lloose_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_lloose_2018", "hist_lloose_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "looseljet_pt", "looseljet_eta");
    auto hist_lnotag_2018 = df1.Filter("sample_year == \"2018\"").Histo2D({"hist_lnotag_2018", "hist_lnotag_2018;jet_pt;jet_eta", 10, pt_bins, 4, eta_bins}, "notagljet_pt", "notagljet_eta");

    TFile* f = new TFile("hist_beff.root", "RECREATE");

    std::cout << "Processing 2016preVFP hists..." << std::endl;

    TH2D hist_btotal_2016preVFP;
    hist_btight_2016preVFP.GetValue().Copy(hist_btotal_2016preVFP);
    hist_btotal_2016preVFP.Add(&(hist_bloose_2016preVFP.GetValue()));
    hist_btotal_2016preVFP.Add(&(hist_bnotag_2016preVFP.GetValue()));

    TH2D hist_ctotal_2016preVFP;
    hist_ctight_2016preVFP.GetValue().Copy(hist_ctotal_2016preVFP);
    hist_ctotal_2016preVFP.Add(&(hist_cloose_2016preVFP.GetValue()));
    hist_ctotal_2016preVFP.Add(&(hist_cnotag_2016preVFP.GetValue()));

    TH2D hist_ltotal_2016preVFP;
    hist_ltight_2016preVFP.GetValue().Copy(hist_ltotal_2016preVFP);
    hist_ltotal_2016preVFP.Add(&(hist_lloose_2016preVFP.GetValue()));
    hist_ltotal_2016preVFP.Add(&(hist_lnotag_2016preVFP.GetValue()));

    hist_btight_2016preVFP->Divide(&hist_btotal_2016preVFP);
    hist_bloose_2016preVFP->Divide(&hist_btotal_2016preVFP);
    hist_bnotag_2016preVFP->Divide(&hist_btotal_2016preVFP);

    hist_ctight_2016preVFP->Divide(&hist_ctotal_2016preVFP);
    hist_cloose_2016preVFP->Divide(&hist_ctotal_2016preVFP);
    hist_cnotag_2016preVFP->Divide(&hist_ctotal_2016preVFP);

    hist_ltight_2016preVFP->Divide(&hist_ltotal_2016preVFP);
    hist_lloose_2016preVFP->Divide(&hist_ltotal_2016preVFP);
    hist_lnotag_2016preVFP->Divide(&hist_ltotal_2016preVFP);

    std::cout << "Writing 2016preVFP hists..." << std::endl;

    hist_btight_2016preVFP->Write();
    hist_bloose_2016preVFP->Write();
    hist_bnotag_2016preVFP->Write();

    hist_ctight_2016preVFP->Write();
    hist_cloose_2016preVFP->Write();
    hist_cnotag_2016preVFP->Write();

    hist_ltight_2016preVFP->Write();
    hist_lloose_2016preVFP->Write();
    hist_lnotag_2016preVFP->Write();

    std::cout << "2016preVFP hists written..." << std::endl;

    std::cout << "Processing 2016postVFP hists..." << std::endl;

   
    TH2D hist_btotal_2016postVFP;
    hist_btight_2016postVFP.GetValue().Copy(hist_btotal_2016postVFP);
    hist_btotal_2016postVFP.Add(&(hist_bloose_2016postVFP.GetValue()));
    hist_btotal_2016postVFP.Add(&(hist_bnotag_2016postVFP.GetValue()));

    TH2D hist_ctotal_2016postVFP;
    hist_ctight_2016postVFP.GetValue().Copy(hist_ctotal_2016postVFP);
    hist_ctotal_2016postVFP.Add(&(hist_cloose_2016postVFP.GetValue()));
    hist_ctotal_2016postVFP.Add(&(hist_cnotag_2016postVFP.GetValue()));

    TH2D hist_ltotal_2016postVFP;
    hist_ltight_2016postVFP.GetValue().Copy(hist_ltotal_2016postVFP);
    hist_ltotal_2016postVFP.Add(&(hist_lloose_2016postVFP.GetValue()));
    hist_ltotal_2016postVFP.Add(&(hist_lnotag_2016postVFP.GetValue()));

    hist_btight_2016postVFP->Divide(&hist_btotal_2016postVFP);
    hist_bloose_2016postVFP->Divide(&hist_btotal_2016postVFP);
    hist_bnotag_2016postVFP->Divide(&hist_btotal_2016postVFP);

    hist_ctight_2016postVFP->Divide(&hist_ctotal_2016postVFP);
    hist_cloose_2016postVFP->Divide(&hist_ctotal_2016postVFP);
    hist_cnotag_2016postVFP->Divide(&hist_ctotal_2016postVFP);

    hist_ltight_2016postVFP->Divide(&hist_ltotal_2016postVFP);
    hist_lloose_2016postVFP->Divide(&hist_ltotal_2016postVFP);
    hist_lnotag_2016postVFP->Divide(&hist_ltotal_2016postVFP);

    std::cout << "Writing 2016postVFP hists..." << std::endl;

    hist_btight_2016postVFP->Write();
    hist_bloose_2016postVFP->Write();
    hist_bnotag_2016postVFP->Write();

    hist_ctight_2016postVFP->Write();
    hist_cloose_2016postVFP->Write();
    hist_cnotag_2016postVFP->Write();

    hist_ltight_2016postVFP->Write();
    hist_lloose_2016postVFP->Write();
    hist_lnotag_2016postVFP->Write();

    std::cout << "2016postVFP hists written..." << std::endl;

    std::cout << "Processing 2017 hists..." << std::endl;

   
    TH2D hist_btotal_2017;
    hist_btight_2017.GetValue().Copy(hist_btotal_2017);
    hist_btotal_2017.Add(&(hist_bloose_2017.GetValue()));
    hist_btotal_2017.Add(&(hist_bnotag_2017.GetValue()));

    TH2D hist_ctotal_2017;
    hist_ctight_2017.GetValue().Copy(hist_ctotal_2017);
    hist_ctotal_2017.Add(&(hist_cloose_2017.GetValue()));
    hist_ctotal_2017.Add(&(hist_cnotag_2017.GetValue()));

    TH2D hist_ltotal_2017;
    hist_ltight_2017.GetValue().Copy(hist_ltotal_2017);
    hist_ltotal_2017.Add(&(hist_lloose_2017.GetValue()));
    hist_ltotal_2017.Add(&(hist_lnotag_2017.GetValue()));

    hist_btight_2017->Divide(&hist_btotal_2017);
    hist_bloose_2017->Divide(&hist_btotal_2017);
    hist_bnotag_2017->Divide(&hist_btotal_2017);

    hist_ctight_2017->Divide(&hist_ctotal_2017);
    hist_cloose_2017->Divide(&hist_ctotal_2017);
    hist_cnotag_2017->Divide(&hist_ctotal_2017);

    hist_ltight_2017->Divide(&hist_ltotal_2017);
    hist_lloose_2017->Divide(&hist_ltotal_2017);
    hist_lnotag_2017->Divide(&hist_ltotal_2017);

    std::cout << "Writing 2017 hists..." << std::endl;

    hist_btight_2017->Write();
    hist_bloose_2017->Write();
    hist_bnotag_2017->Write();

    hist_ctight_2017->Write();
    hist_cloose_2017->Write();
    hist_cnotag_2017->Write();

    hist_ltight_2017->Write();
    hist_lloose_2017->Write();
    hist_lnotag_2017->Write();

    std::cout << "2017 hists written..." << std::endl;

    std::cout << "Processing 2018 hists..." << std::endl;

    
    TH2D hist_btotal_2018;
    hist_btight_2018.GetValue().Copy(hist_btotal_2018);
    hist_btotal_2018.Add(&(hist_bloose_2018.GetValue()));
    hist_btotal_2018.Add(&(hist_bnotag_2018.GetValue()));

    TH2D hist_ctotal_2018;
    hist_ctight_2018.GetValue().Copy(hist_ctotal_2018);
    hist_ctotal_2018.Add(&(hist_cloose_2018.GetValue()));
    hist_ctotal_2018.Add(&(hist_cnotag_2018.GetValue()));

    TH2D hist_ltotal_2018;
    hist_ltight_2018.GetValue().Copy(hist_ltotal_2018);
    hist_ltotal_2018.Add(&(hist_lloose_2018.GetValue()));
    hist_ltotal_2018.Add(&(hist_lnotag_2018.GetValue()));

    hist_btight_2018->Divide(&hist_btotal_2018);
    hist_bloose_2018->Divide(&hist_btotal_2018);
    hist_bnotag_2018->Divide(&hist_btotal_2018);

    hist_ctight_2018->Divide(&hist_ctotal_2018);
    hist_cloose_2018->Divide(&hist_ctotal_2018);
    hist_cnotag_2018->Divide(&hist_ctotal_2018);

    hist_ltight_2018->Divide(&hist_ltotal_2018);
    hist_lloose_2018->Divide(&hist_ltotal_2018);
    hist_lnotag_2018->Divide(&hist_ltotal_2018);

    std::cout << "Writing 2018 hists..." << std::endl;

    hist_btight_2018->Write();
    hist_bloose_2018->Write();
    hist_bnotag_2018->Write();

    hist_ctight_2018->Write();
    hist_cloose_2018->Write();
    hist_cnotag_2018->Write();

    hist_ltight_2018->Write();
    hist_lloose_2018->Write();
    hist_lnotag_2018->Write();
    
    std::cout << "2018 hists written..." << std::endl;

    f->Close();
    return 0;
}