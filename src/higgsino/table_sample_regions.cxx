///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>
#include <math.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLegend.h" // Controls error level reporting
#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/event_scan.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/hist1d.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "core/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  // sample can be search, ttbar, zll, qcd
  string sample_name = "search";
  // year can be 2016, 2017, 2018 or 2016,2017,2018
  string year_string = "2016";
  bool unblind = false;
  bool single_thread = false;
  bool no_signal = false;
  bool no_t5hh = true;
  bool do_htobb = false;
  bool no_mc = false;
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  // Options: split_mc_in_detail,
  // Options: veto_events_with_list,
  // Options: make_event_list,make_event_list_for_excess
  // Options: process_event_list_data,process_event_list_mc,process_event_list_signal
  // Options: sigscan_points,no_weight
  // Options: TChiHH_HToAll
  string string_options = "";
}

string regSearch(string inString, string const & regExPattern) {
    std::smatch match;
    std::regex pattern(regExPattern);
    std::regex_search (inString, match, pattern);
    return match.str();
}

string getPartialDatasetName(string filename) {
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string fullDatasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  //cout<<datasetName.substr(0,datasetName.find_last_of("Tune"))<<endl;;
  string datasetName = fullDatasetName;
  if (datasetName.find("Run") != string::npos && datasetName.find("RunII") == string::npos) datasetName = "Data";
  if (datasetName.find("SMS-TChiHH") != string::npos) {
    if (datasetName.find("2D") != string::npos) datasetName = "TChiHH_HToBB_HToBB_2D_Tune";
    else datasetName = "TChiHH_HToBB_HToBB_Tune";
  }
  if (datasetName.find("Tune") != string::npos) datasetName = datasetName.substr(0, datasetName.find("Tune")+4);
  if (datasetName.find("_13TeV") != string::npos) datasetName = datasetName.substr(0, datasetName.find("_13TeV"));
  return datasetName;
}

pair<int, int> getSignalMassValues(string filename) {
  //pair<int, int> nlsp_lsp_mass;
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string datasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  if (datasetName.find("SMS-TChiHH") == string::npos) return {-1,-1};
  string NLSPMassString = regSearch(datasetName, "mChi-.*?_").substr(5);
  string LSPMassString = regSearch(datasetName, "mLSP-.*?_").substr(5);
  //cout<<NLSPMassString<<" "<<LSPMassString<<endl;
  return {stoi(NLSPMassString),stoi(LSPMassString)};
}

const NamedFunc eventNumberVeto("eventNumberVeto", [](const Baby &b) -> NamedFunc::ScalarType{
  vector<map<Long64_t, set<tuple<string, int, int, int> > > > * eventVetoData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int> > > > *> (b.EventVetoData());
  bool found = false;
  if ((*eventVetoData)[0].count(b.event())==1) {
    for (auto const & it : (*eventVetoData)[0][b.event()]) {
      if (get<1>(it) == abs(b.SampleType()) && get<3>(it) == b.lumiblock() && get<2>(it) == b.run()) {
        if (getPartialDatasetName(*b.FileNames().begin()) == get<0>(it)) {
          //cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
          found = true;
          break;
        }
      }
    }
  }

  //int index = 0;
  //for (auto const & name : b.FileNames()) {
  //  if(index!=0) cout<<index<<" "<<name<<" "<<b.SampleType()<<" "<<b.run()<<" "<<b.lumiblock()<<" "<<b.event()<<endl;
  //  index = index+1;
  //}
  return !found;
});

string getLuminosityString(string const & local_year_string) {
  set<int> years;
  HigUtilities::parseYears(local_year_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();
  return total_luminosity_string;
}

// Sets folders and filenames
// nanoAodFolder: /net/cms29/cms29r0/pico/NanoAODv5/
// production: higgsino_humboldt
// dataType: mc/data/signal
// fileTag: (tt, single_t, zjets, wjets, qcd, other, all), (all), (200_1, 600_1, 950_1, ...)
// sample: search/ttbar/zll/qcd
// year_string: 2016/2017/2018/run2
void addProcess(string const & processName, Process::Type type, int color, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & dataType, string const & fileTag, string const & local_sample_name, string const & local_year_string, 
  vector<shared_ptr<Process> > & procs) {
  set<int> years;
  //years_a = {2016};
  HigUtilities::parseYears(local_year_string, years);

  // Set base tags
  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["ttonly"]     = set<string>({"*TTJets_*Lept*","*_TTGJets*.root"});
  mctags["ttz"]     = set<string>({"*_TTZ*.root"});
  mctags["ttw"]     = set<string>({"*_TTW*.root"});
  mctags["tth"]     = set<string>({"*ttHTobb*.root"});
  mctags["tttt"]     = set<string>({"*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  //mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
  mctags["zjets"]   = set<string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = set<string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = set<string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  set<string> fileNames;
  if (dataType == "data") fileNames = {"*.root"};
  else if (dataType == "signal") {
    vector<string> massPoints;
    HigUtilities::stringToVectorString(fileTag, massPoints, "_");
    // Special case for angeles production
    if (production == "higgsino_angeles" && massPoints[1] == "0") massPoints[1] = "1";
    if (fileTag=="1D") fileNames = {"*TChiHH*HToBB_Tune*.root"};
    else if (fileTag=="2D") fileNames = {"*TChiHH*HToBB_2D_Tune*.root"};
    else fileNames = {"*TChiHH_mChi-"+massPoints[0]+"_mLSP-"+massPoints[1]+"_*.root"};
  }
  else if (dataType=="t5hh") {
    if (fileTag=="1D") fileNames = {"*T5qqqqZH*.root"};
    else if (fileTag=="2D") fileNames = {"*T5qqqqZH_HToBB*.root"};
    else fileNames = set<string>({fileTag});
  }
  else fileNames = mctags[fileTag];

  // Set folders
  torch::OrderedDict<string, string> folderDict;
  string mcProductionFolder = nanoAodFolder+"/"+production;
  string dataProductionFolder = nanoAodFolder+"/"+production;
  string signalProductionFolder = nanoAodFolder+"/"+production;
  string t5hhProductionFolder = nanoAodFolder+"/"+production;
  folderDict.insert("mc_production_folder", mcProductionFolder);
  // preselect:
  //   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
  //   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
  //folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_higloose/");
  // higloose: 
  //   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
  //   met>150 && nvlep==0
  folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_higloose/");
  //folderDict.insert("search_mc_skim_folder", "mc/skim_met150/");
  // higlep1T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5)) &&
  //   nlep==1 && 
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) 
  folderDict.insert("ttbar_mc_skim_folder", "mc/merged_higmc_higlep1T/");
  // higlep2T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nlep==2 && 
  //   @ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) // pass_2l_trig30
  folderDict.insert("zll_mc_skim_folder", "mc/merged_higmc_higlep2T/");
  // higqcd with met150:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nvlep==0 && ntk==0 && low_dphi_met &&
  //   met>150  // Since applied to met150 skim
  folderDict.insert("qcd_mc_skim_folder", "mc/merged_higmc_higqcd/");

  folderDict.insert("data_production_folder", dataProductionFolder);
  folderDict.insert("search_data_skim_folder","data/merged_higdata_higloose/");
  //folderDict.insert("search_data_skim_folder","data/skim_met150/");
  folderDict.insert("ttbar_data_skim_folder","data/merged_higdata_higlep1T/");
  folderDict.insert("zll_data_skim_folder","data/merged_higdata_higlep2T/");
  folderDict.insert("qcd_data_skim_folder","data/merged_higdata_higqcd/");

  folderDict.insert("signal_production_folder", signalProductionFolder);
  if (HigUtilities::is_in_string_options(string_options, "TChiHH_HToAll")) folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_HToAll_fastSimJmeCorrection/merged_higmc_higloose/");
  else folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higloose/");
  //folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_2D/skim_met150/");
  folderDict.insert("ttbar_signal_skim_folder", "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep1T/");
  folderDict.insert("zll_signal_skim_folder", "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higlep2T/");
  folderDict.insert("qcd_signal_skim_folder", "SMS-TChiHH_2D_fastSimJmeCorrection/merged_higmc_higqcd/");

  folderDict.insert("t5hh_production_folder", t5hhProductionFolder);
  if (dataType=="t5hh" && fileTag=="2D") {
    folderDict.insert("search_t5hh_skim_folder", "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higloose/");
    folderDict.insert("ttbar_t5hh_skim_folder", "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higlep1T/");
    folderDict.insert("zll_t5hh_skim_folder", "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higlep2T/");
    folderDict.insert("qcd_t5hh_skim_folder", "SMS-T5qqqqZH_fastSimJmeCorrection/merged_higmc_higqcd/");
  }
  else {
    folderDict.insert("search_t5hh_skim_folder", "SMS-T5qqqqZH_FullSim/merged_higmc_higloose/");
    folderDict.insert("ttbar_t5hh_skim_folder", "SMS-T5qqqqZH_FullSim/merged_higmc_higlep1T/");
    folderDict.insert("zll_t5hh_skim_folder", "SMS-T5qqqqZH_FullSim/merged_higmc_higlep2T/");
    folderDict.insert("qcd_t5hh_skim_folder", "SMS-T5qqqqZH_FullSim/merged_higmc_higqcd/");
  }

  // Set paths
  set<string> pathNames;
  pathNames = attach_folder(folderDict[dataType+"_production_folder"], years, folderDict[local_sample_name+"_"+dataType+"_skim_folder"], fileNames);

  if (dataType == "data") procs.push_back(Process::MakeShared<Baby_pico>(processName, type, color, pathNames, additionalCut));
  else procs.push_back(Process::MakeShared<Baby_pico>(processName, type, color, pathNames,"stitch"&&additionalCut)); // mc, signal
}

void addAllMcProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & local_sample_name, string const & local_year_string, 
  vector<shared_ptr<Process> > & procs) {
  Palette colors("txt/colors.txt", "default");
  addProcess("t#bar{t}+X (#tau_{had}>0) "+processName_postfix, Process::Type::background, colors("tt_htau"), "ntrutauh>0"&&additionalCut,
    nanoAodFolder, production, "mc", "tt", local_sample_name, local_year_string,
    procs);
  addProcess("t#bar{t}+X (#tau_{had}=0) "+processName_postfix, Process::Type::background, colors("tt_1l"), "ntrutauh==0"&&additionalCut,
    nanoAodFolder, production, "mc", "tt", local_sample_name, local_year_string,
    procs);
  addProcess("Z+Jets "+processName_postfix, Process::Type::background, kOrange+1, additionalCut,
    nanoAodFolder, production, "mc", "zjets", local_sample_name, local_year_string,
    procs);
  addProcess("W+Jets "+processName_postfix, Process::Type::background, kGreen+1, additionalCut,
    nanoAodFolder, production, "mc", "wjets", local_sample_name, local_year_string,
    procs);
  addProcess("Single t "+processName_postfix, Process::Type::background, colors("single_t"), additionalCut,
    nanoAodFolder, production, "mc", "single_t", local_sample_name, local_year_string,
    procs);
  addProcess("QCD "+processName_postfix, Process::Type::background, colors("other"), additionalCut,
    nanoAodFolder, production, "mc", "qcd", local_sample_name, local_year_string,
    procs);
  addProcess("Other "+processName_postfix, Process::Type::background, kGray+2, additionalCut,
    nanoAodFolder, production, "mc", "other", local_sample_name, local_year_string,
    procs);
}

void addAllMcDetailProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & local_sample_name, string const & local_year_string, 
  vector<shared_ptr<Process> > & procs) {
  Palette colors("txt/colors.txt", "default");
  addProcess("t#bar{t} "+processName_postfix, Process::Type::background, colors("tt_htau"), additionalCut,
    nanoAodFolder, production, "mc", "ttonly", local_sample_name, local_year_string,
    procs);
  addProcess("t#bar{t}W "+processName_postfix, Process::Type::background, colors("tt_htau"), additionalCut,
    nanoAodFolder, production, "mc", "ttw", local_sample_name, local_year_string,
    procs);
  addProcess("t#bar{t}Z "+processName_postfix, Process::Type::background, colors("tt_htau"), additionalCut,
    nanoAodFolder, production, "mc", "ttz", local_sample_name, local_year_string,
    procs);
  addProcess("t#bar{t}H "+processName_postfix, Process::Type::background, colors("tt_htau"), additionalCut,
    nanoAodFolder, production, "mc", "tth", local_sample_name, local_year_string,
    procs);
  addProcess("tt#bar{t}#bar{t} "+processName_postfix, Process::Type::background, colors("tt_htau"), additionalCut,
    nanoAodFolder, production, "mc", "tttt", local_sample_name, local_year_string,
    procs);
  addProcess("Z+Jets "+processName_postfix, Process::Type::background, kOrange+1, additionalCut,
    nanoAodFolder, production, "mc", "zjets", local_sample_name, local_year_string,
    procs);
  addProcess("W+Jets "+processName_postfix, Process::Type::background, kGreen+1, additionalCut,
    nanoAodFolder, production, "mc", "wjets", local_sample_name, local_year_string,
    procs);
  addProcess("Single t "+processName_postfix, Process::Type::background, colors("single_t"), additionalCut,
    nanoAodFolder, production, "mc", "single_t", local_sample_name, local_year_string,
    procs);
  addProcess("QCD "+processName_postfix, Process::Type::background, colors("other"), additionalCut,
    nanoAodFolder, production, "mc", "qcd", local_sample_name, local_year_string,
    procs);
  addProcess("Other "+processName_postfix, Process::Type::background, kGray+2, additionalCut,
    nanoAodFolder, production, "mc", "other", local_sample_name, local_year_string,
    procs);
}

void addMultipleSignalProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & local_sample_name, string const & local_year_string, 
  vector<shared_ptr<Process> > & procs) 
{
  //special case for checking bright point on significance scan
  if (HigUtilities::is_in_string_options(string_options, "sigscan_points")) {
    vector<string> sigm = {"250", "275", "300"};
    vector<string> sigmlsp = {"100", "100", "100"};
    for (unsigned isig(0); isig<sigm.size(); isig++){
      addProcess("TChiHH("+sigm[isig]+","+sigmlsp[isig]+") "+processName_postfix, Process::Type::signal, 1, additionalCut,
        nanoAodFolder, production, "signal", sigm[isig]+"_"+sigmlsp[isig], local_sample_name, local_year_string,
        procs);
    }
  }
  else {
    vector<string> sigm = {"200", "500", "950"};
    for (unsigned isig(0); isig<sigm.size(); isig++){
      addProcess("TChiHH("+sigm[isig]+",0) "+processName_postfix, Process::Type::signal, 1, additionalCut,
        nanoAodFolder, production, "signal", sigm[isig]+"_0", local_sample_name, local_year_string,
        procs);
    }
    if (!no_t5hh) {
      vector<string> sigm_t5hh = {"*mGluino-1000_mChi-950_mLSP-1_*.root","*mGluino-1600_mChi-1550_mLSP-1_*.root","*mGluino-2000_mChi-1950_mLSP-1_*.root"};
      vector<string> model_names = {"T5HH(1000,0)","T5HH(1600,0)","T5HH(2000,0)"};
      //vector<string> sigm_t5hh = {"*mGluino-1200_mChi-1150_mLSP-400_*.root","*mGluino-1600_mChi-1550_mLSP-1_*.root","*mGluino-2000_mChi-1950_mLSP-1_*.root"};
      //vector<string> model_names = {"T5HH(1200_400)","T5HH(1600_0)","T5HH(2000_0)"};
      for (unsigned isig(0); isig<sigm_t5hh.size(); isig++){
        addProcess(model_names[isig]+processName_postfix, Process::Type::signal, 1, additionalCut,
            nanoAodFolder, production, "t5hh", sigm_t5hh[isig], local_sample_name, local_year_string,
            procs);
      }
    }
  }
}

void addDataProcess(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & sample, string const & local_year_string, 
  vector<shared_ptr<Process> > & procs) 
{
    //NamedFunc lepton_triggers = "(HLT_IsoMu24 || HLT_IsoMu27 || HLT_Mu50 || HLT_Ele27_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT)";
    //NamedFunc met_triggers = "(HLT_PFMET110_PFMHT110_IDTight || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight || HLT_PFMET120_PFMHT120_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFMET120_PFMHT120_IDTight_PFHT60 || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)";
    NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
    NamedFunc met_triggers = Higfuncs::met_trigger;;

    NamedFunc triggers_data = "1";
    if (sample == "zll") triggers_data = lepton_triggers;
    else if (sample == "ttbar") triggers_data = lepton_triggers || met_triggers;
    else if (sample == "qcd") triggers_data = met_triggers;
    else  triggers_data = met_triggers;

    addProcess("Data"+processName_postfix, Process::Type::data, kBlack, triggers_data && additionalCut,
      nanoAodFolder, production, "data", "", sample, local_year_string,
      procs);
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);
  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  // dataType: mc/data/signal
  // fileTag: (tt, single_t, zjets, wjets, qcd, other, all), (all), (200_1, 600_1, 950_1, ...)
  // sample: search/ttbar/zll/qcd
  // year_string: 2016/2017/2018/run2
  string production_a = "higgsino_klamath"; 
  string nanoAodFolder_a;
  nanoAodFolder_a = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7";
  string sample_a = sample_name;
  string year_string_a = year_string;

  string total_luminosity_string = HigUtilities::getLuminosityString(year_string_a);

  //    Define processes, including intersections
  //--------------------------------------------------
  NamedFunc base_filters = "1";
  if (sample_a == "search") base_filters = Higfuncs::final_pass_filters;
  else if (sample_a == "ttbar") base_filters = Higfuncs::final_ttbar_pass_filters;
  else if (sample_a == "zll") base_filters = Higfuncs::final_zll_pass_filters;
  else if (sample_a == "qcd") base_filters = Higfuncs::final_qcd_pass_filters;

  NamedFunc search_resolved_cuts = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved_cuts = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         && Higfuncs::lead_signal_lepton_pt>30;
  NamedFunc zll_resolved_cuts =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==2&&njet>=4&&njet<=5&&met<50&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  NamedFunc qcd_resolved_cuts =
                         "met/mht<2 && met/met_calo<2&&"
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";

  NamedFunc resolved_cuts = "1";
  if (sample_a == "ttbar") resolved_cuts = ttbar_resolved_cuts;
  else if (sample_a == "zll") resolved_cuts = zll_resolved_cuts;
  else if (sample_a == "qcd") resolved_cuts = qcd_resolved_cuts;
  else resolved_cuts = search_resolved_cuts;

  vector<shared_ptr<Process> > procs;
  if (!no_mc) {
    bool split_mc_in_detail = HigUtilities::is_in_string_options(string_options, "split_mc_in_detail");
    if (split_mc_in_detail) {
      addAllMcDetailProcesses("", base_filters&&resolved_cuts,
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs);
    } else {
      addAllMcProcesses("", base_filters&&resolved_cuts,
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs);
    }
  }
  if (!no_signal) {addMultipleSignalProcesses("", base_filters&&resolved_cuts,
    nanoAodFolder_a, production_a, 
    sample_a, year_string_a,
    procs);
  }
  if (unblind) {
    addDataProcess("", base_filters&&resolved_cuts,
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs);
  }

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  // For event list
  vector<shared_ptr<Process> > procs_all;
  bool make_event_list = HigUtilities::is_in_string_options(string_options, "make_event_list");
  if (make_event_list) {
    if (years.count(2016)==1) {
      if (!no_mc) {addProcess("MC_2016", Process::Type::background, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "mc", "all", sample_a, "2016",
          procs_all);
      }
      if (!no_signal) {addProcess("TChiHH1D_2016", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "signal", "1D", sample_a, "2016",
          procs_all);
        addProcess("TChiHH2D_2016", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "signal", "2D", sample_a, "2016",
          procs_all);
        addProcess("T5HH1D_2016", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "t5hh", "1D", sample_a, "2016",
          procs_all);
        addProcess("T5HH2D_2016", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "t5hh", "2D", sample_a, "2016",
          procs_all);
      }
      if (unblind) {
        addDataProcess("_2016", base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, 
          sample_a, "2016",
          procs_all);
      }
    }
    if (years.count(2017)==1) {
      if (!no_mc) {addProcess("MC_2017", Process::Type::background, 1, base_filters&&resolved_cuts,
        nanoAodFolder_a, production_a, "mc", "all", sample_a, "2017",
        procs_all);
      }
      if (!no_signal) {addProcess("TChiHH1D_2017", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "signal", "1D", sample_a, "2017",
          procs_all);
        addProcess("TChiHH2D_2017", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "signal", "2D", sample_a, "2017",
          procs_all);
        addProcess("T5HH1D_2017", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "t5hh", "1D", sample_a, "2017",
          procs_all);
        addProcess("T5HH2D_2017", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "t5hh", "2D", sample_a, "2017",
          procs_all);
      }
      if (unblind) {
        addDataProcess("_2017", base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, 
          sample_a, "2017",
          procs_all);
      }
    }
    if (years.count(2018)==1) {
      if (!no_mc) {addProcess("MC_2018", Process::Type::background, 1, base_filters&&resolved_cuts,
        nanoAodFolder_a, production_a, "mc", "all", sample_a, "2018",
        procs_all);
      }
      if (!no_signal) {addProcess("TChiHH1D_2018", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "signal", "1D", sample_a, "2018",
          procs_all);
        addProcess("TChiHH2D_2018", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "signal", "2D", sample_a, "2018",
          procs_all);
        addProcess("T5HH1D_2018", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "t5hh", "1D", sample_a, "2018",
          procs_all);
        addProcess("T5HH2D_2018", Process::Type::signal, 1, base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, "t5hh", "2D", sample_a, "2018",
          procs_all);
      }
      if (unblind) {
        addDataProcess("_2018", base_filters&&resolved_cuts,
          nanoAodFolder_a, production_a, 
          sample_a, "2018",
          procs_all);
      }
    }
  }

  NamedFunc weight = Higfuncs::final_weight;
  if (do_htobb) 
    weight = Higfuncs::final_weight*Higfuncs::htobb;
  if (HigUtilities::is_in_string_options(string_options, "sigscan_points"))
    weight *= HigUtilities::w_CNToN1N2;
  if (HigUtilities::is_in_string_options(string_options, "no_weight")) 
    weight = "1"; 

  //    Useful binning definitions
  //------------------------------------------
  bool bin_met = true;
  bool bin_drmax = true;

  vector<TString> vc_drmax;
  vc_drmax.push_back("hig_cand_drmax[0]<=1.1");
  vc_drmax.push_back("hig_cand_drmax[0]>1.1");
  if (!bin_drmax) vc_drmax = {"hig_cand_drmax[0]>0"};
  vector<string> res_abcd;
  vector<TString> vc_met;

  string c_2b, c_3b, c_4b, hig, sbd;
  if (sample_a=="qcd") {
    c_2b = "nbm==0";
    c_3b = "nbm==1";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig};

    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=400");
    vc_met.push_back("met>400");
    if (!bin_met) vc_met = {"met>150"};
  } else if (sample_a=="zll") {
    c_2b = "nbm==0";
    c_3b = "nbm==1";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig};

    vc_met.push_back("ll_pt[0]>0   && ll_pt[0]<=75");
    vc_met.push_back("ll_pt[0]>75  && ll_pt[0]<=150");
    vc_met.push_back("ll_pt[0]>150 && ll_pt[0]<=200");
    vc_met.push_back("ll_pt[0]>200 && ll_pt[0]<=300");
    vc_met.push_back("ll_pt[0]>300");
    if (!bin_met) vc_met = {"ll_pt[0]>0"};
  } else if (sample_a=="ttbar") {
    c_2b = "(nbt==2&&nbm==2)";
    c_3b = "(nbt>=2&&nbm==3&&nbl==3)";
    c_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig, c_4b+"&&"+sbd, c_4b+"&&"+hig};

    vc_met.push_back("met>0   && met<=75");
    vc_met.push_back("met>75  && met<=150");
    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300");
    if (!bin_met) vc_met = {"met>0"};
  } else { // search
    c_2b = "(nbt==2&&nbm==2)";
    c_3b = "(nbt>=2&&nbm==3&&nbl==3)";
    c_4b = "(nbt>=2&&nbm>=3&&nbl>=4)";
    hig = "hig_cand_am[0]>100 && hig_cand_am[0]<=140";
    sbd = "!(hig_cand_am[0]>100 && hig_cand_am[0]<=140)";
    res_abcd = {c_2b+"&&"+sbd, c_2b+"&&"+hig, c_3b+"&&"+sbd, c_3b+"&&"+hig, c_4b+"&&"+sbd, c_4b+"&&"+hig};

    vc_met.push_back("met>150 && met<=200");
    vc_met.push_back("met>200 && met<=300");
    vc_met.push_back("met>300 && met<=400");
    vc_met.push_back("met>400");
    if (!bin_met) vc_met = {"met>150"};
  }


  //        Define and fill tables
  //------------------------------------
  { // To be able to destruct pm so that files are written to system before going to next step.
    PlotMaker pm;
    vector<TableRow> tablerows;
    unsigned nbins(0);
    for (auto &imet: vc_met) {
      for (auto &idrmax: vc_drmax) {
        // Label
        tablerows.push_back(TableRow("$"+CodeToLatex((imet+"&&"+idrmax).Data())+"$"));
        for (auto &iabcd: res_abcd) {
          // Cut for row
          NamedFunc res_reg_ = imet+"&&"+idrmax+"&&"+iabcd;
          // Add row to table
          tablerows.push_back(TableRow("$"+CodeToLatex(iabcd)+"$", res_reg_, 0,0, weight));
          nbins++;
        }
      }
    }
    TString tabname = "table_"; tabname += "resolved";

    bool print_uncertainty = true;
    pm.Push<Table>("FixName:table_"+sample_a+"_regions", tablerows, procs, 0, 1, 0, 1, 0, print_uncertainty).LuminosityTag(total_luminosity_string);

    // signal region
    NamedFunc signalRegion = "hig_cand_am[0]>100 && hig_cand_am[0]<140 && nbm>=3";

    // Load event veto data
    bool veto_events_with_list = HigUtilities::is_in_string_options(string_options, "veto_events_with_list");
    vector<map<Long64_t, set<tuple<string, int, int, int> > > > * eventVetoData = new vector<map<Long64_t, set<tuple<string, int, int, int> > > >(1, map<Long64_t, set<tuple<string, int, int, int> > > ());
    if (veto_events_with_list) {
      set<string> vetoEventListNames = {"processed_resolved_list_CR_SCAN_MC.txt", "processed_resolved_list_SR_SCAN_MC.txt",
                                        "processed_resolved_list_CR_SCAN_TChiHH.txt", "processed_resolved_list_SR_SCAN_TChiHH.txt"};
      //eventVetoData[0][event] = {[sampleType, year, run, lumiblock]}
      for (string const & vetoEventListName : vetoEventListNames) {
        ifstream vetoEventList(vetoEventListName);
        string line;
        while (getline(vetoEventList, line)) {
          if(line[0]=='#') continue;
          vector<string> lineSplit;
          HigUtilities::stringToVectorString(line, lineSplit, ",");
          string const & sampleType = lineSplit[0];
          int year = stoi(lineSplit[1]);
          int run = stoi(lineSplit[2]);
          int lumiblock = stoi(lineSplit[3]);
          Long64_t event = stoll(lineSplit[4]);
          // Insert data
          if ((*eventVetoData)[0].count(event) == 0) {
            (*eventVetoData)[0][event] = {{sampleType, year, run, lumiblock}};
          } else {
            if ((*eventVetoData)[0][event].count({sampleType, year, run, lumiblock}) != 0) cout<<"Duplicate: "<<sampleType<<" "<<year<<" "<<run<<" "<<lumiblock<<" "<<event<<endl;
            (*eventVetoData)[0][event].insert({sampleType, year, run, lumiblock});
          }
        }
        vetoEventList.close();
      }
    }

    pm.SetEventVetoData(static_cast<void * >(eventVetoData));

    pm.min_print_ = true;
    pm.multithreaded_ = !single_thread;

    if (make_event_list) {
      pm.Push<EventScan>("resolved_list_SR", base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event", "met"}, procs_all, 10, true);
      pm.Push<EventScan>("resolved_list_CR", base_filters&&search_resolved_cuts&&!signalRegion,vector<NamedFunc>{"SampleType", "run","lumiblock","event", "met"}, procs_all, 10, true);
      //pm.Push<EventScan>("resolved_list_SR_met300", "met>300"&&base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event", "met"}, procs_all, 10, true);
      //pm.Push<EventScan>("resolved_list_CR_met300", "met>300"&&base_filters&&search_resolved_cuts&&!signalRegion,vector<NamedFunc>{"SampleType", "run","lumiblock","event", "met"}, procs_all, 10, true);
      ////pm.Push<EventScan>("resolved_list_SR_met300", base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event", "met"}, procs_all, 10, true);
      ////pm.Push<EventScan>("resolved_list_CR_met300", base_filters&&search_resolved_cuts&&!signalRegion,vector<NamedFunc>{"SampleType", "run","lumiblock","event", "met"}, procs_all, 10, true);
      //pm.Push<EventScan>("resolved_list_SR_met300", "met>300"&&base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event", "met"}, procs_all, 10, true);
      //pm.Push<EventScan>("resolved_list_CR_met300", "met>300"&&base_filters&&search_resolved_cuts&&!signalRegion,vector<NamedFunc>{"SampleType", "run","lumiblock","event", "met"}, procs_all, 10, true);
    }

    bool make_event_list_for_excess = HigUtilities::is_in_string_options(string_options, "make_event_list_for_excess");
    if (make_event_list_for_excess) {
      pm.Push<EventScan>("SR_met200_lowdrmax", "met>200&&hig_cand_drmax[0]<=1.1"&&base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event", "met", "hig_cand_drmax[0]","hig_cand_am[0]", "hig_cand_dm[0]"}, procs_all, 10, true);
    }

    pm.MakePlots(1.);
  } // destruct pm so that files are written to system before going to next step.

  bool process_event_list_data = HigUtilities::is_in_string_options(string_options, "process_event_list_data");
  bool process_event_list_mc = HigUtilities::is_in_string_options(string_options, "process_event_list_mc");
  bool process_event_list_signal = HigUtilities::is_in_string_options(string_options, "process_event_list_signal");
  set<string> eventListNames = {};
  if (process_event_list_data) {
    //eventListNames.insert("resolved_list_CR_met300_SCAN_Data_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_Data_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_Data_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_Data_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_Data_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_Data_2018.txt");
    eventListNames.insert("resolved_list_CR_SCAN_Data_2016.txt");eventListNames.insert("resolved_list_SR_SCAN_Data_2016.txt");
    eventListNames.insert("resolved_list_CR_SCAN_Data_2017.txt");eventListNames.insert("resolved_list_SR_SCAN_Data_2017.txt");
    eventListNames.insert("resolved_list_CR_SCAN_Data_2018.txt");eventListNames.insert("resolved_list_SR_SCAN_Data_2018.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_Data_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_Data_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_Data_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_Data_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_Data_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_Data_2018.txt");
  }
  if (process_event_list_mc) {
    //eventListNames.insert("resolved_list_CR_met300_SCAN_MC_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_MC_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_MC_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_MC_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_MC_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_MC_2018.txt");
    eventListNames.insert("resolved_list_CR_SCAN_MC_2016.txt");eventListNames.insert("resolved_list_SR_SCAN_MC_2016.txt");
    eventListNames.insert("resolved_list_CR_SCAN_MC_2017.txt");eventListNames.insert("resolved_list_SR_SCAN_MC_2017.txt");
    eventListNames.insert("resolved_list_CR_SCAN_MC_2018.txt");eventListNames.insert("resolved_list_SR_SCAN_MC_2018.txt");
  }
  if (process_event_list_signal) {
    //eventListNames.insert("resolved_list_CR_met300_SCAN_TChiHH1D_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_TChiHH1D_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_TChiHH1D_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_TChiHH1D_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_TChiHH1D_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_TChiHH1D_2018.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_TChiHH2D_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_TChiHH2D_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_TChiHH2D_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_TChiHH2D_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_TChiHH2D_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_TChiHH2D_2018.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_T5HH1D_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_T5HH1D_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_T5HH1D_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_T5HH1D_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_T5HH1D_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_T5HH1D_2018.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_T5HH2D_2016.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_T5HH2D_2016.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_T5HH2D_2017.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_T5HH2D_2017.txt");
    //eventListNames.insert("resolved_list_CR_met300_SCAN_T5HH2D_2018.txt");eventListNames.insert("resolved_list_SR_met300_SCAN_T5HH2D_2018.txt");
    eventListNames.insert("resolved_list_CR_SCAN_TChiHH1D_2016.txt");eventListNames.insert("resolved_list_SR_SCAN_TChiHH1D_2016.txt");
    eventListNames.insert("resolved_list_CR_SCAN_TChiHH1D_2017.txt");eventListNames.insert("resolved_list_SR_SCAN_TChiHH1D_2017.txt");
    eventListNames.insert("resolved_list_CR_SCAN_TChiHH1D_2018.txt");eventListNames.insert("resolved_list_SR_SCAN_TChiHH1D_2018.txt");
    eventListNames.insert("resolved_list_CR_SCAN_TChiHH2D_2016.txt");eventListNames.insert("resolved_list_SR_SCAN_TChiHH2D_2016.txt");
    eventListNames.insert("resolved_list_CR_SCAN_TChiHH2D_2017.txt");eventListNames.insert("resolved_list_SR_SCAN_TChiHH2D_2017.txt");
    eventListNames.insert("resolved_list_CR_SCAN_TChiHH2D_2018.txt");eventListNames.insert("resolved_list_SR_SCAN_TChiHH2D_2018.txt");
    eventListNames.insert("resolved_list_CR_SCAN_T5HH1D_2016.txt");eventListNames.insert("resolved_list_SR_SCAN_T5HH1D_2016.txt");
    eventListNames.insert("resolved_list_CR_SCAN_T5HH1D_2017.txt");eventListNames.insert("resolved_list_SR_SCAN_T5HH1D_2017.txt");
    eventListNames.insert("resolved_list_CR_SCAN_T5HH1D_2018.txt");eventListNames.insert("resolved_list_SR_SCAN_T5HH1D_2018.txt");
    eventListNames.insert("resolved_list_CR_SCAN_T5HH2D_2016.txt");eventListNames.insert("resolved_list_SR_SCAN_T5HH2D_2016.txt");
    eventListNames.insert("resolved_list_CR_SCAN_T5HH2D_2017.txt");eventListNames.insert("resolved_list_SR_SCAN_T5HH2D_2017.txt");
    eventListNames.insert("resolved_list_CR_SCAN_T5HH2D_2018.txt");eventListNames.insert("resolved_list_SR_SCAN_T5HH2D_2018.txt");
  }
  for (string const & eventListName : eventListNames) {
    ifstream eventList(eventListName);
    bool isSignal = false;
    if (eventListName.find("TChiHH") != string::npos) isSignal = true;
    if (eventListName.find("T5HH") != string::npos) isSignal = true;
    string processedEventListName = "processed_"+eventListName;
    cout<<"Processing "<<eventListName<<" to "<<processedEventListName<<endl;
    ofstream processedEventList(processedEventListName);
    if (isSignal) processedEventList<<"# SampleName, Year, RunNumber, LumiBlockNumber, EventNumber, MET, NLSP mass, LSP mass"<<endl;
    else processedEventList<<"# SampleName, Year, RunNumber, LumiBlockNumber, EventNumber, MET"<<endl;
    string line;
    while (getline(eventList, line)) {
      if (Contains(line, "Instance")) continue;
      vector<string> lineSplit;
      HigUtilities::stringToVectorString(line, lineSplit, " ");
      // Expect "Row Instance SampleType run lumiblock event Filename" sequence
      int sampleType = stoi(lineSplit[2]); string run = lineSplit[3]; string lumiblock = lineSplit[4]; string event = lineSplit[5]; string met_string = lineSplit[6]; string filename = lineSplit[7];
      int met = round(stof(met_string));
      string datasetName = getPartialDatasetName(filename);
      if (isSignal) {
        pair<int, int> nlsp_lsp_mass = getSignalMassValues(filename);
        processedEventList<<datasetName<<", "<<abs(sampleType)<<", "<<run<<", "<<lumiblock<<", "<<event<<", "<<met<<", "<<nlsp_lsp_mass.first<<", "<<nlsp_lsp_mass.second<<endl;
      } else processedEventList<<datasetName<<", "<<abs(sampleType)<<", "<<run<<", "<<lumiblock<<", "<<event<<", "<<met<<endl;
    }
    eventList.close();
    processedEventList.close();
  };

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 0},
      {"sample", required_argument, 0, 's'},
      {"year", required_argument, 0, 'y'},
      {"unblind", no_argument, 0, 'u'},
      {"no_signal", no_argument, 0, 0},
      {"htobb", no_argument, 0, 0},
      {"t5hh", no_argument, 0, 0},
      {"no_mc", no_argument, 0, 0},
      {"string_options", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    // put : for required_argument
    opt = getopt_long(argc, argv, "s:y:uo:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      sample_name = optarg;
      break;
    case 'y':
      year_string = optarg;
      break;
    case 'u':
      unblind = true;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "single_thread"){
        single_thread = true;
      } else if (optname == "no_signal") {
        no_signal = true;
      } else if (optname == "no_mc") {
        no_mc = true;
      } else if (optname == "t5hh") {
        no_t5hh = false;
      } else if (optname == "htobb") {
        do_htobb = true;
      } else {
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
