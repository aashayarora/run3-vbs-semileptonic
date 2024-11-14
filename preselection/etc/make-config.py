from glob import glob
import json
import re
import argparse
import uproot

class Config:
    def __init__(self, sample_category : str, samples: str):
        self.sample_category = sample_category
        self.samples = sorted(glob(samples))
        self.config = {"samples": {}}

        self.process_samples(xsecs)
        self.write_config(f"{sample_category}.json")

    def write_config(self, output_file):
        with open(output_file, "w") as f:
            json.dump(self.config, f, indent=4)

    @staticmethod
    def extract_sample_year(sample):
        if "UL16" in sample and "APV" in sample:
            return "2016preVFP"
        elif "UL16" in sample and not "APV" in sample:
            return "2016postVFP"
        elif "UL17" in sample or "UL2017" in sample:
            return "2017"
        elif "UL18" in sample or "UL2018" in sample:
            return "2018"
        elif any(run in sample for run in ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F"]) and "HIPM" in sample:
            return "2016preVFP"
        elif any(run in sample for run in ["Run2016F", "Run2016G", "Run2016H"]) and not "HIPM" in sample:
            return "2016postVFP"
        else:
            raise ValueError(f"Error: year not found for {sample}")

    @staticmethod
    def get_xsec_weight(xsecs, sample):
        if sample in xsecs:
            return xsecs[sample]
        else:
            raise ValueError(f"xsec not found for {sample}")

    @staticmethod
    def get_lumi(year):
        lumi = {
            "2016preVFP": 19.52,
            "2016postVFP": 16.81,
            "2017": 41.529,
            "2018": 59.7
        }
        if year in lumi:
            return lumi[year]
        else:
            raise ValueError(f"lumi not found for {year}")

    @staticmethod
    def extract_mc_sample_type(sample_name):
        sample_type_mapping = {
            "DY": "DY",
            "TTTo": "ttbar",
            "TT": "ttx",
            "tt": "ttx",
            "ST": "ST",
            "WJets": "WJets",
            "EWK": "EWK"
        }
        for key, value in sample_type_mapping.items():
            if key in sample_name:
                return value
        return "Other"

    @staticmethod
    def get_sample_name(sample):
        if "data" in sample:
            return re.search(r"/([^/]+)_Run20", sample).group(1)
        else:
            return re.search(r"/([^/]+)_TuneCP5", sample).group(1)

    def process_samples(self, xsecs):
        for sample in self.samples:
            try:
                sample_name = self.get_sample_name(sample)
                sample_year = self.extract_sample_year(sample)
                xsec = self.get_xsec_weight(xsecs, f"{sample_name},{sample_year}") if self.sample_category != "data" else 1.0
                num_events = 0
                if self.sample_category != "data":
                    files = glob(f"{sample}/output*.root")
                    for file in files:
                        with uproot.open(file) as upf:
                            num_events += sum(upf["Runs"]["genEventSumw"].array())
                else:
                    num_events = 1.0
                self.config["samples"].update(
                    {
                        f"{sample_name}_{sample_year}": {
                            "trees": ["Events"],
                            "files": [f"{sample}/output*.root"],
                            "metadata": {
                                "sample_category": self.sample_category,
                                "sample_year": sample_year,
                                "sample_type": self.extract_mc_sample_type(sample_name) if self.sample_category != "data" else "SingleMuon" if "SingleMuon" in sample_name else "SingleElectron",
                                "xsec": xsec,
                                "lumi": self.get_lumi(sample_year) if self.sample_category != "data" else 1.0,
                                "nevents": num_events
                            }
                        }
                    }
                )
            except Exception as e:
                print(f"Error in {sample}: {e}")

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--category", help="categories: bkg, sig or data", required=True)
    args = argparser.parse_args()

    if not args.category or not any(cat in args.category for cat in ["bkg", "sig", "data"]):
        raise ValueError("Please provide a valid category")

    with open("xsecs.json", "r") as f_xsecs:
        xsecs = json.load(f_xsecs)

    skim_paths = {
        "bkg": "/ceph/cms/store/user/aaarora/VBS_1lep_skims/bkg_1lep_4ak4_or_1ak8_2ak4_v1/*",
        "data": "/ceph/cms/store/user/aaarora/VBS_1lep_skims/data_1lep_4ak4_or_1ak8_2ak4_v1/*",
        "sig": "/ceph/cms/store/user/aaarora/VBS_1lep_skims/sig_1lep_4ak4_or_1ak8_2ak4_v1/*/*/NANOAODSIM/*/*/skimmed/*.root",
        "sig_private": "/ceph/cms/store/user/jguiang/VBSVHSkim/sig_1lep_4ak4_or_1ak8_2ak4_v1/*Inclusive*/*.root"
    }

    config = Config(args.category, skim_paths[args.category])