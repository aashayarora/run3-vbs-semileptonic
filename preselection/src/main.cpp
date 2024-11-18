#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

#include "weights.h"
#include "selections.h"
#include "corrections.h"
#include "utils.h"

#include "argparser.hpp"

struct MyArgs : public argparse::Args {
    std::string &spec = kwarg("i,input", "spec.json path");
    bool &cutflow = flag("cutflow", "print cutflow");
    bool &JEC = flag("jec", "JEC");
    bool &JER = flag("jer", "JER");
    bool &JMS = flag("jms", "JMS");
    bool &JMR = flag("jmr", "JMR");
    bool &METUnclustered = flag("met", "MET unclustered");
    int &nthreads = kwarg("n,nthreads", "number of threads").set_default(1);
    std::string &output = kwarg("o,output", "output root file").set_default("");
    std::string &variation = kwarg("var", "variation").set_default("nominal");
    std::string &JEC_type = kwarg("jec_type", "JEC type").set_default("");
    std::string &cut = kwarg("cut", "cut on final snapshot").set_default("passCut9");
};

void runDataAnalysis(RNode df, MyArgs args, std::string output_file) {
    auto df_ = defineCorrectedCols(df);
    auto df1 = ObjectSelections(df_);
    auto df_metcorr = METPhiCorrections(cset_met_2016preVFP, cset_met_2016postVFP, cset_met_2017, cset_met_2018, df1);
    auto df_event = EventSelections(df_metcorr);
    auto df_weights = applyDataWeights(df_event);
    
    std::vector<std::string> cuts = {"passCut1", "passCut2", "passCut3", "passCut4", "passCut5", "passCut6", "passCut7", "passCut8", "passCut9", "passCut8_cr", "passCut9_cr"};
    auto cutflow = Cutflow(df_weights, cuts);

    df_weights = df_weights.Filter(args.cut);
    saveSnapshot(df_weights, std::string(output_file), true);

    if (args.cutflow) {
        cutflow.Print("cutflow/" + output_file + "_cutflow.txt");
    }
}

void runMCAnalysis(RNode df, MyArgs args, std::string output_file) {
    // corrections
    auto df_ = defineCorrectedCols(df);
    auto df1 = ObjectSelections(df_);
    // apply pre preselection corrections
    df1 = METPhiCorrections(cset_met_2016preVFP, cset_met_2016postVFP, cset_met_2017, cset_met_2018, df1);
    if (args.JER) {
        df1 = JetEnergyResolution(cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, cset_jer_smear, df1, args.variation);
    } else if (args.METUnclustered) {
        df1 = METUnclusteredCorrections(df1, args.variation);
    } else if (args.JEC) {
        df1 = JetEnergyCorrection(cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, df1, args.JEC_type, args.variation);
    }
    
    auto df_event = EventSelections(df1);

    if (args.JMS) {
        df_event = JMS_Corrections(cset_jms, df_event, args.variation);
    } else if (args.JMR) {
        df_event = JMR_Corrections(cset_jmr, df_event, args.variation);
    }
    auto df_weights = applyMCWeights(df_event);
    
    std::vector<std::string> cuts = {"passCut1", "passCut2", "passCut3", "passCut4", "passCut5", "passCut6", "passCut7", "passCut8", "passCut9", "passCut8_cr", "passCut9_cr"};
    auto cutflow = Cutflow(df_weights, cuts);
    
    df_weights = df_weights.Filter(args.cut);
    saveSnapshot(df_weights, std::string(output_file));

    if (args.cutflow){
        cutflow.Print("cutflow/" + output_file + "_cutflow.txt");
    }
}

int main(int argc, char** argv) {
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;

    ROOT::EnableImplicitMT(args.nthreads);
    ROOT::RDataFrame df_ = ROOT::RDF::Experimental::FromSpec(input_spec);
    ROOT::RDF::Experimental::AddProgressBar(df_);

    // define metadata
    auto df = defineMetadata(df_);
    // run analysis
    if (input_spec.find("data") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "data";
        }
        runDataAnalysis(df, args, output_file);
    }
    else {
        if (output_file.empty()) {
            if (input_spec.find("bkg") != std::string::npos) {
                output_file = "bkg";
            }
            else if (input_spec.find("sig") != std::string::npos) {
                output_file = "sig";
            }
            else {
                std::cerr << "Incorrect spec name, file must contain sig, bkg or data" << std::endl;
            }
        }
        runMCAnalysis(df, args, output_file);
    }
    return 0;
}