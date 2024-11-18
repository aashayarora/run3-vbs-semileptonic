#include "selections.h"

RNode METFilters(RNode df) {
    return df.Define("passesFlags", "Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter");
}

RNode Triggers(RNode df) {
    return df.Define("passesTriggers", "");
}

RNode ElectronSelections(RNode df) {
    return df.Define("Electron_SC_absEta", "Electron_eta + Electron_deltaEtaSC")
            .Define("vetoElectrons", 
                "Electron_pt > 7 &&"
                "abs(Electron_SC_absEta) < 2.5 && "
                "abs(Electron_dxy) <= 0.05 && "
                "abs(Electron_dz) < 0.1 && "
                "abs(Electron_sip3d) < 8 && "
                "Electron_miniPFRelIso_all < 0.4 && "
                "Electron_lostHits <= 1 && "
                "Electron_mvaFall17V2noIso_WPL == true")
            .Define("nVetoElectrons", "nElectron == 0 ? 0 : Sum(vetoElectrons)")
            .Define("electron_jetBTagDeepFlav", mediumDFBtagWP, {"sample_year"})
            .Define("tightElectrons", "vetoElectrons &&" 
                "Electron_mvaTTHUL > 0.9 && "
                "Electron_pt > 35 && "
                "Electron_hoe < 0.1 && "
                "Electron_eInvMinusPInv > -0.04 && "
                "((abs(Electron_SC_absEta) <= 1.479 && Electron_sieie < 0.011) || ((abs(Electron_SC_absEta) > 1.479 && abs(Electron_SC_absEta) <= 2.5) && Electron_sieie <= 0.030)) && "
                "Electron_convVeto == true && "
                "Electron_tightCharge == 2 && "
                "Electron_lostHits == 0 && "
                "Electron_jetBTagDeepFlav < electron_jetBTagDeepFlav")
            .Define("nTightElectrons", "nElectron == 0 ? 0 : Sum(tightElectrons)")
            .Define("electron_pt", "Electron_pt[tightElectrons]")
            .Define("electron_eta", "Electron_eta[tightElectrons]")
            .Define("electron_phi", "Electron_phi[tightElectrons]")
            .Define("electron_mass", "Electron_mass[tightElectrons]")
            .Define("isElectron", "nVetoElectrons == 1 && nTightElectrons == 1");
}

RNode MuonSelections(RNode df) {
    return df.Define("vetoMuons", 
                "Muon_pt > 5 && "
                "abs(Muon_eta) < 2.4 && "
                "abs(Muon_dxy) < 0.05 && "
                "abs(Muon_dz) < 0.1 && "
                "abs(Muon_sip3d) < 8 && "
                "Muon_miniPFRelIso_all <= 0.4 && "
                "Muon_looseId == 1")
            .Define("nVetoMuons", "nMuon == 0 ? 0 : Sum(vetoMuons)")
            .Define("muon_jetBTagDeepFlav", mediumDFBtagWP, {"sample_year"})
            .Define("tightMuons", "vetoMuons && "
                "Muon_pt > 35 && "
                "Muon_mediumId && "
                "Muon_mvaTTHUL > 0.85 && "
                "Muon_jetBTagDeepFlav < muon_jetBTagDeepFlav")
            .Define("nTightMuons", "nMuon == 0 ? 0 : Sum(tightMuons)")
            .Define("muon_pt", "Muon_pt[tightMuons]")
            .Define("muon_eta", "Muon_eta[tightMuons]")
            .Define("muon_phi", "Muon_phi[tightMuons]")
            .Define("muon_mass", "Muon_mass[tightMuons]")
            .Define("isMuon", "nVetoMuons == 1 && nTightMuons == 1");
}

RNode LeptonSelections(RNode df) {
    auto df_el = electronSelections(df);
    auto df_mu = muonSelections(df_el);
    return df_mu.Define("lepton_pt", "isElectron ? GElectron_pt[0] : GMuon_pt[0]")
            .Define("lepton_eta", "isElectron ? GElectron_eta[0] : GMuon_eta[0]")
            .Define("lepton_phi", "isElectron ? GElectron_phi[0] : GMuon_phi[0]")
            .Define("lepton_mass", "isElectron ? GElectron_mass[0] : GMuon_mass[0]");
}

RNode HbbSelections(RNode df) {
    return df.Define("HLepDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "GLepton_eta", "GLepton_phi"})
            .Define("HScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
            .Define("HCandidateJets", 
                "CorrFatJet_pt > 250 && "
                "HLepDeltaR >= 0.8 && "
                "CorrFatJet_mass > 50 && "
                "abs(FatJet_eta) <= 2.5 && "
                "FatJet_msoftdrop > 40 && "
                "FatJet_jetId > 0")
            .Define("HighestHScoreIdx", "HScore.size() != 0 ? ArgMax(HScore[HCandidateJets]) : 999.0")
            .Define("hbb_score", "HighestHScoreIdx != 999.0 ? HScore[HCandidateJets][HighestHScoreIdx] : -1.0")
            .Define("hbb_pt", "CorrFatJet_pt[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_eta", "FatJet_eta[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_phi", "FatJet_phi[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_mass", "FatJet_particleNet_mass[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_msoftdrop", "FatJet_msoftdrop[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_tau2", "FatJet_tau2[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_tau1", "FatJet_tau1[HCandidateJets][HighestHScoreIdx]")
            .Define("hbb_tau21", "hbb_tau2 / hbb_tau1");
}

RNode VqqSelections(RNode df) {
    return df.Define("WLepDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "GLepton_eta", "GLepton_phi"})
            .Define("WHDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "GHiggs_eta", "GHiggs_phi"})
            .Define("WZCandidateJets", 
                "CorrFatJet_pt > 250 && "
                "CorrFatJet_mass > 50 && "
                "WLepDeltaR >= 0.8 && "
                "WHDeltaR >= 0.8 && "
                "abs(FatJet_eta) <= 2.5 && "
                "FatJet_msoftdrop > 40 && "
                "FatJet_jetId > 0")
            .Define("WScore", "(FatJet_particleNetMD_Xqq + FatJet_particleNetMD_Xcc) / (FatJet_particleNetMD_Xqq + FatJet_particleNetMD_Xcc + FatJet_particleNetMD_QCD)")
            .Define("ZScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
            .Define("HighestWjetScoreIdx", "WScore.size() != 0 ? ArgMax(WScore[WZCandidateJets]) : -1")
            .Define("HighestZjetScoreIdx", "ZScore.size() != 0 ? ArgMax(ZScore[WZCandidateJets]) : -1")
            .Define("vzqq_score", "HighestZjetScoreIdx != -1 ? ZScore[WZCandidateJets][HighestZjetScoreIdx] : -1")
            .Define("vqq_score", "HighestWjetScoreIdx != -1 ? WScore[WZCandidateJets][HighestWjetScoreIdx] : -1")
            .Define("vqq_pt", "CorrFatJet_pt[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_eta", "FatJet_eta[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_phi", "FatJet_phi[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_mass", "FatJet_particleNet_mass[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_msoftdrop", "FatJet_msoftdrop[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_tau2", "FatJet_tau2[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_tau1", "FatJet_tau1[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("vqq_tau21", "vqq_tau2 / vqq_tau1");
}

RNode AK4Selections(RNode df) {
    return df.Define("AK4LepDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "GLepton_eta", "GLepton_phi"})
            .Define("AK4HDeltaR", VfDeltaR, { "Jet_eta", "Jet_phi", "GHiggs_eta", "GHiggs_phi"})
            .Define("AK4WDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "GW_eta", "GW_phi"})
            .Define("ak4tightBjetScore", tightDFBtagWP, {"sample_year"})
            .Define("ak4looseBjetScore", looseDFBtagWP, {"sample_year"})
            .Define("goodJets", "CorrJet_pt >= 20 && "
                "abs(Jet_eta) < 2.5 && "
                "AK4LepDeltaR >= 0.4 && "
                "AK4HDeltaR >= 0.8 && "
                "AK4WDeltaR >= 0.8 && "
                "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
            .Define("ak4FromBJet", "goodJets && Jet_btagDeepFlavB > ak4tightBjetScore")
            .Define("goodLooseBJets", "goodJets && Jet_btagDeepFlavB > ak4looseBjetScore")
            .Define("sortedBJets", "Argsort(-Jet_btagDeepFlavB[goodLooseBJets])")
            .Define("GBJet_pt", "Take(CorrJet_pt[goodLooseBJets], sortedBJets)")
            .Define("GBJet_eta", "Take(Jet_eta[goodLooseBJets], sortedBJets)")
            .Define("GBJet_phi", "Take(Jet_phi[goodLooseBJets], sortedBJets)")
            .Define("GBJet_mass", "Take(CorrJet_mass[goodLooseBJets], sortedBJets)")
            .Define("GBJet_score", "Take(Jet_btagDeepFlavB[goodLooseBJets], sortedBJets)")
            .Define("bjet1pt", "GBJet_pt.size() > 0 ? GBJet_pt[0] : -999")
            .Define("bjet1eta", "GBJet_pt.size() > 0 ? GBJet_eta[0] : -999")
            .Define("bjet1phi", "GBJet_pt.size() > 0 ? GBJet_phi[0] : -999")
            .Define("bjet1score", "GBJet_pt.size() > 0 ? GBJet_score[0] : -999")
            .Define("bjet2pt", "GBJet_pt.size() > 1 ? GBJet_pt[1] : -999")
            .Define("bjet2eta", "GBJet_pt.size() > 1 ? GBJet_eta[1] : -999")
            .Define("bjet2phi", "GBJet_pt.size() > 1 ? GBJet_phi[1] : -999")
            .Define("bjet2score", "GBJet_pt.size() > 1 ? GBJet_score[1] : -999")
            .Define("Mlb", VfInvariantMass, {"GBJet_pt", "GBJet_eta", "GBJet_phi", "GBJet_mass", "GLepton_pt", "GLepton_eta", "GLepton_phi", "GLepton_mass"})
            .Define("MinMlbJetIdx", "Mlb.size() != 0 ? ArgMin(Mlb) : -1")
            .Define("Mlbminloose", "Mlb.size() != 0 ? Mlb[MinMlbJetIdx] : 1000");
}

RNode VBSJetsSelections(RNode df) {
    return df.Define("goodVBSJets", "CorrJet_pt >= 30 && "
                "abs(Jet_eta) <= 4.7 && "
                "AK4LepDeltaR >= 0.4 && "
                "AK4HDeltaR >= 0.8 && "
                "AK4WDeltaR >= 0.8 && "
                "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
            .Define("VBSJets_pt", "CorrJet_pt[goodVBSJets]")
            .Define("VBSJets_eta", "Jet_eta[goodVBSJets]")
            .Define("VBSJets_phi", "Jet_phi[goodVBSJets]")
            .Define("VBSJets_mass", "CorrJet_mass[goodVBSJets]")
            .Define("VBSjetidxs", VBS_MaxE, {"VBSJets_pt", "VBSJets_eta", "VBSJets_phi", "VBSJets_mass"})
            .Define("vbs1_pt", "VBSJets_pt[VBSjetidxs[0]]")
            .Define("vbs1_eta", "VBSJets_eta[VBSjetidxs[0]]")
            .Define("vbs1_phi", "VBSJets_phi[VBSjetidxs[0]]")
            .Define("vbs1_mass", "VBSJets_mass[VBSjetidxs[0]]")
            .Define("vbs2_pt", "VBSJets_pt[VBSjetidxs[1]]")
            .Define("vbs2_eta", "VBSJets_eta[VBSjetidxs[1]]")
            .Define("vbs2_phi", "VBSJets_phi[VBSjetidxs[1]]")
            .Define("vbs2_mass", "VBSJets_mass[VBSjetidxs[1]]")
            .Define("vbs_ptjj", "VBSjet1pt + VBSjet2pt")
            .Define("vbs_detajj", "abs(VBSjet1eta - VBSjet2eta)")
            .Define("vbs_mjj", fInvariantMass, {"VBSjet1pt", "VBSjet1eta", "VBSjet1phi", "VBSjet1mass", "VBSjet2pt", "VBSjet2eta", "VBSjet2phi", "VBSjet2mass"})
            .Define("met", "CorrMET_pt")
            .Define("ST", "lepton_pt + met + hbb_pt + vqq_pt");
}

RNode ObjectSelections(RNode df) {
    auto df_flags = METFilters(df);
    auto df_trigger = TriggerSelections(df_flags);
    auto df_lepton = LeptonSelections(df_trigger);
    auto df_higgs = HbbSelections(df_lepton);
    auto df_w = VqqSelections(df_higgs);
    auto df_ak4 = AK4Selections(df_w);
    auto df_vbs = VBSJetsSelections(df_ak4);
    return df_vbs; 
}

RNode EventSelections(RNode df) {
     return df_vbs.Define("passCut1", "passesFlags && passesTriggers")
        .Define("passCut2", "passCut1 && ((nVetoMuons == 1 && nTightMuons == 1 && nVetoElectrons == 0 && nTightElectrons == 0) || "
            "(nVetoMuons == 0 && nTightMuons == 0 && nVetoElectrons == 1 && nTightElectrons == 1)) && "
            "(GLepton_pt > 40)")
        .Define("passCut3", "passCut2 && HighestHScore > 0")
        .Define("passCut4", "passCut3 && HighestWjetScore > 0")
        .Define("passCut5", "passCut4 && Sum(ak4FromBJet) == 0")
        .Define("passCut6", "passCut5 && Sum(goodVBSJets) >= 2")
        .Define("passCut7", "passCut6 && ST > 1000")
        .Define("passCut8", "passCut7 && HighestHScore > 0.5")
        .Define("passCut9", "passCut8 && HighestWjetScore > 0.7")
        .Define("passCut8_cr", "passCut7 && HighestHScore < 0.95")
        .Define("passCut9_cr", "passCut8_cr && HighestWjetScore < 0.7");
}