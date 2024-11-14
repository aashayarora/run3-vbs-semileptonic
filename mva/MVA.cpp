#include "MVA.h"

struct MyArgs : public argparse::Args {
    std::string &input = kwarg("i,input", "input path");
    std::string &output = kwarg("o,output", "output path").set_default("");
    std::string &cut = kwarg("cut", "cut on final snapshot").set_default("passCut9");
    int &nthreads = kwarg("n,nthreads", "number of threads").set_default(1);
    bool &old = flag("old", "use old MVA");
};

int main(int argc, char** argv){
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_file = args.input;

    ROOT::EnableImplicitMT(args.nthreads);
    ROOT::RDataFrame df("Events", input_file);
    ROOT::RDF::Experimental::AddProgressBar(df);

    std::string output_file = args.output.empty() ? input_file.substr(0, input_file.find_last_of('.')) + "_MVA.root" : args.output;

    if (args.old) {
        TMVA::Experimental::RReader reader_AB("weights/BDT/AB/TMVAClassification_BDT.weights.xml");
        TMVA::Experimental::RReader reader_BA("weights/BDT/BA/TMVAClassification_BDT.weights.xml");

        TMVA::Experimental::RSofieReader dnn("weights/DNN/old_model.pt", {{1, 7}});

        auto predict = [&](unsigned long long event, float VBSjet1pt, float VBSjet1eta, float VBSjet1phi, float VBSjet2pt, float VBSjet2eta, float VBSjet2phi, float VBSMjj, float VBSdetajj){
            if(event % 2){
                auto score_AB = reader_AB.Compute({VBSjet1pt, VBSjet1eta, VBSjet1phi, VBSjet2pt, VBSjet2eta, VBSjet2phi, VBSMjj, VBSdetajj});
                return score_AB[0];
            }
            else{
                auto score_BA = reader_BA.Compute({VBSjet1pt, VBSjet1eta, VBSjet1phi, VBSjet2pt, VBSjet2eta, VBSjet2phi, VBSMjj, VBSdetajj});
                return score_BA[0];
            }
        };

        auto df1 = df.Filter(args.cut)
                .Define("VBSBDTscore", predict, {"event", "VBSjet1pt", "VBSjet1eta", "VBSjet1phi", "VBSjet2pt","VBSjet2eta", "VBSjet2phi", "VBSMjj", "VBSdetajj"})
                .Define("DNN_Hbbmass", "(float) (Hbbmass - 50) / 200")
                .Define("DNN_Wjetmass", "(float) Wjetmass / 200")
                .Define("DNN_HbbPt", "(float) log(HbbPt)")
                .Define("DNN_WjetPt", "(float) log(WjetPt)")
                .Define("DNN_leptonpt", "(float) log(leptonpt)")
                .Define("DNN_Mlbminloose", "(float) (Mlbminloose - 0) / 1000")
                .Define("DNN_MET", "(float) log(MET)")
                .Define("abcdnet_output", Compute<7, float>(dnn), {"DNN_Hbbmass", "DNN_HbbPt", "DNN_Wjetmass", "DNN_WjetPt",  "DNN_leptonpt", "DNN_Mlbminloose", "DNN_MET"})
                .Define("abcdnet_score", "abcdnet_output[0]");

        df1.Snapshot("Events", output_file);

        return 0;

    } else {
        TMVA::Experimental::RBDT bdt("VBSBDT", "weights/BDT/BDT_Weights.root");
        TMVA::Experimental::RSofieReader dnn("weights/DNN/model.pt", {{1, 7}});

        auto df1 = df.Filter(args.cut)
                .Define("VBSBDTOutput", Compute<8, float>(bdt), {"VBSjet1pt", "VBSjet1eta", "VBSjet1phi", "VBSjet2pt","VBSjet2eta", "VBSjet2phi", "VBSMjj", "VBSdetajj"})
                .Define("VBSBDTscore", "VBSBDTOutput[0]")
                .Define("DNN_Hbbmass", "(float) (Hbbmass - 50) / 200")
                .Define("DNN_Wjetmass", "(float) Wjetmass / 200")
                .Define("DNN_HbbPt", "(float) log(HbbPt)")
                .Define("DNN_WjetPt", "(float) log(WjetPt)")
                .Define("DNN_leptonpt", "(float) log(leptonpt)")
                .Define("DNN_Mlbminloose", "(float) (Mlbminloose - 0) / 1000")
                .Define("DNN_MET", "(float) log(MET)")
                .Define("abcdnet_output", Compute<7, float>(dnn), {"DNN_Hbbmass", "DNN_HbbPt", "DNN_Wjetmass", "DNN_WjetPt",  "DNN_leptonpt", "DNN_Mlbminloose", "DNN_MET"})
                .Define("abcdnet_score", "abcdnet_output[0]");

        df1.Snapshot("Events", output_file);

        return 0;
    }
}