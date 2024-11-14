#include "selections.h"

RNode flagSelections(RNode df) {
    auto df_flags = df.Define("passesFlags", "Flag_goodVertices && "
                "(!isData || Flag_globalSuperTightHalo2016Filter) && "
                "Flag_HBHENoiseFilter && "
                "Flag_HBHENoiseIsoFilter && "
                "Flag_EcalDeadCellTriggerPrimitiveFilter && "
                "Flag_BadPFMuonFilter &&"
                "Flag_BadPFMuonDzFilter && "
                "Flag_eeBadScFilter && "
                "(!(is2017 || is2018) || Flag_ecalBadCalibFilter)");
    return df_flags;
}

RNode triggerSelections(RNode df) {
    auto df_triggers = df.Define("passesTriggers", "((!isData) && "
                "((is2018 && (HLT_Ele32_WPTight_Gsf == true || HLT_IsoMu24 == true)) || "
                "(is2017 && (HLT_Ele32_WPTight_Gsf_L1DoubleEG == true || HLT_IsoMu27 == true)) || "
                "(is2016 && (HLT_Ele27_eta2p1_WPTight_Gsf == true || HLT_IsoMu24 == true || HLT_IsoTkMu24 == true)))) || "
                "(isData && "
                "((is2018 && ((sample_type == \"SingleElectron\" && HLT_Ele32_WPTight_Gsf == true) || (sample_type == \"SingleMuon\" && HLT_IsoMu24 == true))) || "
                "(is2017 && ((sample_type == \"SingleElectron\" && HLT_Ele32_WPTight_Gsf_L1DoubleEG == true) || (sample_type == \"SingleMuon\" && HLT_IsoMu27 == true))) || "
                "(is2016 && ((sample_type == \"SingleElectron\" && HLT_Ele27_eta2p1_WPTight_Gsf == true) || (sample_type == \"SingleMuon\" && (HLT_IsoMu24 == true || HLT_IsoTkMu24 == true))))))");
    return df_triggers;
}

RNode electronSelections(RNode df) {
    auto df_el = df.Define("SC_absEta", "Electron_eta + Electron_deltaEtaSC")
            .Define("vetoElectrons", 
                "Electron_pt > 7 &&"
                "abs(SC_absEta) < 2.5 && "
                "abs(Electron_dxy) <= 0.05 && "
                "abs(Electron_dz) < 0.1 && "
                "abs(Electron_sip3d) < 8 && "
                "Electron_miniPFRelIso_all < 0.4 && "
                "Electron_lostHits <= 1 && "
                "Electron_mvaFall17V2noIso_WPL == true")
            .Define("nVetoElectrons", "nElectron == 0 ? 0 : Sum(vetoElectrons)")
            .Define("electronBjetScore", mediumDFBtagWP, {"sample_year"})
            .Define("tightElectrons", "vetoElectrons &&" 
                "Electron_mvaTTHUL > 0.9 && "
                "Electron_pt > 35 && "
                "Electron_hoe < 0.1 && "
                "Electron_eInvMinusPInv > -0.04 && "
                "((abs(SC_absEta) <= 1.479 && Electron_sieie < 0.011) || ((abs(SC_absEta) > 1.479 && abs(SC_absEta) <= 2.5) && Electron_sieie <= 0.030)) && "
                "Electron_convVeto == true && "
                "Electron_tightCharge == 2 && "
                "Electron_lostHits == 0 && "
                "Electron_jetBTagDeepFlav < electronBjetScore")
            .Define("nTightElectrons", "nElectron == 0 ? 0 : Sum(tightElectrons)")
            .Define("GElectron_pt", "Electron_pt[tightElectrons]")
            .Define("GElectron_eta", "Electron_eta[tightElectrons]")
            .Define("GElectron_phi", "Electron_phi[tightElectrons]")
            .Define("GElectron_mass", "Electron_mass[tightElectrons]")
            .Define("isElectron", "nVetoElectrons == 1 && nTightElectrons == 1");
    return df_el;
}

RNode muonSelections(RNode df) {
    auto df_mu = df.Define("vetoMuons", 
                "Muon_pt > 5 && "
                "abs(Muon_eta) < 2.4 && "
                "abs(Muon_dxy) < 0.05 && "
                "abs(Muon_dz) < 0.1 && "
                "abs(Muon_sip3d) < 8 && "
                "Muon_miniPFRelIso_all <= 0.4 && "
                "Muon_looseId == 1")
            .Define("nVetoMuons", "nMuon == 0 ? 0 : Sum(vetoMuons)")
            .Define("muonBjetScore", mediumDFBtagWP, {"sample_year"})
            .Define("tightMuons", "vetoMuons && "
                "Muon_pt > 35 && "
                "Muon_mediumId && "
                "Muon_mvaTTHUL > 0.85 && "
                "Muon_jetBTagDeepFlav < muonBjetScore")
            .Define("nTightMuons", "nMuon == 0 ? 0 : Sum(tightMuons)")
            .Define("GMuon_pt", "Muon_pt[tightMuons]")
            .Define("GMuon_eta", "Muon_eta[tightMuons]")
            .Define("GMuon_phi", "Muon_phi[tightMuons]")
            .Define("GMuon_mass", "Muon_mass[tightMuons]")
            .Define("isMuon", "nVetoMuons == 1 && nTightMuons == 1");
    return df_mu;
}

RNode leptonSelections(RNode df) {
    auto df_el = electronSelections(df);
    auto df_mu = muonSelections(df_el);
    auto df_leps = df_mu.Define("GLepton_pt", "isElectron ? GElectron_pt[0] : GMuon_pt[0]")
            .Define("GLepton_eta", "isElectron ? GElectron_eta[0] : GMuon_eta[0]")
            .Define("GLepton_phi", "isElectron ? GElectron_phi[0] : GMuon_phi[0]")
            .Define("GLepton_mass", "isElectron ? GElectron_mass[0] : GMuon_mass[0]")
            .Define("leptonpt", "GLepton_pt");
    return df_leps;
}

RNode higgsSelections(RNode df) {
    auto df_higgs = df.Define("HLepDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "GLepton_eta", "GLepton_phi"})
            .Define("HScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
            .Define("HCandidateJets", 
                "CorrFatJet_pt > 250 && "
                "HLepDeltaR >= 0.8 && "
                "CorrFatJet_mass > 50 && "
                "abs(FatJet_eta) <= 2.5 && "
                "FatJet_msoftdrop > 40 && "
                "FatJet_jetId > 0")
            .Define("HighestHScoreIdx", "HScore.size() != 0 ? ArgMax(HScore[HCandidateJets]) : 999.0")
            .Define("HighestHScore", "HighestHScoreIdx != 999.0 ? HScore[HCandidateJets][HighestHScoreIdx] : -1.0")
            .Define("GHiggs_pt", "CorrFatJet_pt[HCandidateJets][HighestHScoreIdx]")
            .Define("GHiggs_eta", "FatJet_eta[HCandidateJets][HighestHScoreIdx]")
            .Define("GHiggs_phi", "FatJet_phi[HCandidateJets][HighestHScoreIdx]")
            .Define("GHiggs_mass", "FatJet_particleNet_mass[HCandidateJets][HighestHScoreIdx]")
            .Define("Hbbmass", "FatJet_particleNet_mass[HCandidateJets][HighestHScoreIdx]")
            .Define("Hbbscore", "HighestHScore")
            .Define("HbbPt", "CorrFatJet_pt[HCandidateJets][HighestHScoreIdx]")
            .Define("Htau2", "FatJet_tau2[HCandidateJets][HighestHScoreIdx]")
            .Define("Htau1", "FatJet_tau1[HCandidateJets][HighestHScoreIdx]")
            .Define("Htau21", "Htau2 / Htau1")
            .Define("Hsoftdrop", "FatJet_msoftdrop[HCandidateJets][HighestHScoreIdx]");
    return df_higgs;
}

RNode WZSelections(RNode df) {
    auto df_wz = df.Define("WLepDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "GLepton_eta", "GLepton_phi"})
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
            .Define("HighestWjetScoreIdx", "WScore.size() != 0 ? ArgMax(WScore[WZCandidateJets]) : -1")
            .Define("HighestWjetScore", "HighestWjetScoreIdx != -1 ? WScore[WZCandidateJets][HighestWjetScoreIdx] : -1")
            .Define("ZScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
            .Define("HighestZjetScoreIdx", "ZScore.size() != 0 ? ArgMax(ZScore[WZCandidateJets]) : -1")
            .Define("HighestZjetScore", "HighestZjetScoreIdx != -1 ? ZScore[WZCandidateJets][HighestZjetScoreIdx] : -1")
            .Define("GW_pt", "CorrFatJet_pt[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("GW_eta", "FatJet_eta[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("GW_phi", "FatJet_phi[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("GW_mass", "FatJet_particleNet_mass[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("Wjetmass", "FatJet_particleNet_mass[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("WjetPt", "CorrFatJet_pt[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("Wtau2", "FatJet_tau2[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("Wtau1", "FatJet_tau1[WZCandidateJets][HighestWjetScoreIdx]")
            .Define("Wtau21", "Wtau2 / Wtau1")
            .Define("Wsoftdrop", "FatJet_msoftdrop[WZCandidateJets][HighestWjetScoreIdx]");
    return df_wz;
}

RNode AK4Selections(RNode df) {
    auto df_ak4 = df.Define("AK4LepDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "GLepton_eta", "GLepton_phi"})
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
            .Define("LnTBJets", "goodJets && Jet_btagDeepFlavB > ak4looseBjetScore && Jet_btagDeepFlavB <= ak4tightBjetScore")
            .Define("LnTBJet_pt", "CorrJet_pt[LnTBJets]")
            .Define("LnTBJet_eta", "Jet_eta[LnTBJets]")
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
            // .Define("Ptlb", VfInvariantPt, {"GBJet_pt", "GBJet_eta", "GBJet_phi", "GBJet_mass", "GLepton_pt", "GLepton_eta", "GLepton_phi", "GLepton_mass"})
            // .Define("Philb", VfInvariantPhi, {"GBJet_pt", "GBJet_eta", "GBJet_phi", "GBJet_mass", "GLepton_pt", "GLepton_eta", "GLepton_phi", "GLepton_mass"})
            // .Define("MT", MT, {"Ptlb", "Philb", "CorrMET_pt", "MET_phi"})
            // .Define("MinMT", "MT.size() != 0 ? Min(MT) : 1000")
            .Define("Mlbminloose", "Mlb.size() != 0 ? Mlb[MinMlbJetIdx] : 1000");
    return df_ak4;
}

RNode VBSJetsSelections(RNode df) {
    // auto df_vbs = df.Define("MlbDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "MlbJeteta", "MlbJetphi"})
    auto df_vbs = df.Define("goodVBSJets", "CorrJet_pt >= 30 && "
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
            .Define("puIDJets_pt", "CorrJet_pt[goodVBSJets || goodJets]")
            .Define("puIDJets_eta", "Jet_eta[goodVBSJets || goodJets]")
            .Define("VBSjetidxs", VBS_MaxE, {"VBSJets_pt", "VBSJets_eta", "VBSJets_phi", "VBSJets_mass"})
            .Define("VBSjet1pt", "VBSJets_pt[VBSjetidxs[0]]")
            .Define("VBSjet1eta", "VBSJets_eta[VBSjetidxs[0]]")
            .Define("VBSjet1phi", "VBSJets_phi[VBSjetidxs[0]]")
            .Define("VBSjet1mass", "VBSJets_mass[VBSjetidxs[0]]")
            .Define("VBSjet2pt", "VBSJets_pt[VBSjetidxs[1]]")
            .Define("VBSjet2eta", "VBSJets_eta[VBSjetidxs[1]]")
            .Define("VBSjet2phi", "VBSJets_phi[VBSjetidxs[1]]")
            .Define("VBSjet2mass", "VBSJets_mass[VBSjetidxs[1]]")
            .Define("VBSptjj", "VBSjet1pt + VBSjet2pt")
            .Define("VBSdetajj", "abs(VBSjet1eta - VBSjet2eta)")
            .Define("VBSMjj", fInvariantMass, {"VBSjet1pt", "VBSjet1eta", "VBSjet1phi", "VBSjet1mass", "VBSjet2pt", "VBSjet2eta", "VBSjet2phi", "VBSjet2mass"})
            .Define("MET", "CorrMET_pt")
            .Define("ST", "GLepton_pt + GHiggs_pt + GW_pt + CorrMET_pt");
    return df_vbs;
}

RNode applyPreSelections(RNode df) {
    auto df_flags = flagSelections(df);
    auto df_trigger = triggerSelections(df_flags);
    auto df_lepton = leptonSelections(df_trigger);
    auto df_higgs = higgsSelections(df_lepton);
    auto df_w = WZSelections(df_higgs);
    auto df_ak4 = AK4Selections(df_w);
    auto df_vbs = VBSJetsSelections(df_ak4);
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