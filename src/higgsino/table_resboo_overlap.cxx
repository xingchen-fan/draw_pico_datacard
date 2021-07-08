///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>
#include <regex>

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
#include "higgsino/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  // sample can be search, ttbar, zll, qcd
  string sample_name = "search";
  // year can be 2016, 2017, 2018 or 2016,2017,2018
  string year_string = "2016";
  bool unblind = false;
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
  if (datasetName.find("Run") != string::npos) datasetName = "DATA";
  if (datasetName.find("SMS-TChiHH") != string::npos) datasetName = "TChiHH";
  if (datasetName.find("Tune") != string::npos) datasetName = datasetName.substr(0, datasetName.find("Tune")+4);
  if (datasetName.find("_13TeV") != string::npos) datasetName = datasetName.substr(0, datasetName.find("_13TeV"));
  return datasetName;
}

const NamedFunc eventNumberVeto("eventNumberVeto", [](const Baby &b) -> NamedFunc::ScalarType{
  //b.SampleType();
  //vector<set<tuple<string, int, int, int, Long64_t> > > * eventVetoData = static_cast<vector<set<tuple<string, int, int, int, Long64_t> > > *> (b.EventVetoData());
  //cout<<get<0>((*(*eventVetoData)[0].begin()))<<endl;
  //if ((*eventVetoData)[0].count({})!=0) continue;

  //vector<map<Long64_t, tuple<string, int, int, int> > > * eventVetoData = static_cast<vector<map<Long64_t, tuple<string, int, int, int> > > *> (b.EventVetoData());
  //for (auto const & it : (*eventVetoData)[0]) {
  //  cout<<it.first<<" "<<get<0>(it.second)<<endl;
  //}

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

bool findStringIC(const std::string & strHaystack, const std::string & strNeedle)
{
  auto it = std::search(
    strHaystack.begin(), strHaystack.end(),
    strNeedle.begin(),   strNeedle.end(),
    [](char ch1, char ch2) { return std::toupper(ch1) == std::toupper(ch2); }
  );
  return (it != strHaystack.end() );
}

pair<int, int> getSignalMassValues(string filename) {
  pair<int, int> nlsp_lsp_mass;
  string baseFilename = filename.substr(filename.find_last_of("/\\")+1);
  string datasetName = regSearch(baseFilename, "[A-Z].*|ttHTobb.*");
  if (datasetName.find("SMS-TChiHH") == string::npos) return {-1,-1};
  string NLSPMassString = regSearch(datasetName, "mChi-.*?_").substr(5);
  string LSPMassString = regSearch(datasetName, "mLSP-.*?_").substr(5);
  //cout<<NLSPMassString<<" "<<LSPMassString<<endl;
  return {stoi(NLSPMassString),stoi(LSPMassString)};
}

bool regionCut(const Baby & b, int regionIndex) {
  bool inRegion = false;
  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int, int, int, int> > > > *> (b.EventVetoData());
  // 0: sampleType, 1: year, 2: run, 3: lumiblock, 4: met, 5: nlsp_mass, 6: lsp_mass

  //// Check event nuber
  //if ((*eventNumberData)[regionIndex].count(b.event())==1) {
  //  for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
  //    // Check year and run number
  //    if (get<1>(it) == abs(b.SampleType()) && get<2>(it) == b.run()) {
  //      // Check sample name
  //      if (findStringIC(*b.FileNames().begin(), get<0>(it))) {
  //        // Check nlsp_mass
  //        //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
  //        inRegion = true;
  //        break;
  //      }
  //    }
  //  }
  //}

  // Check event nuber
  if ((*eventNumberData)[regionIndex].count(b.event())==1) {
    for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
      // Check year
      if (get<1>(it) != abs(b.SampleType())) continue;
      // Check run
      if (get<2>(it) != b.run()) continue;
      bool isSignal = ((*b.FileNames().begin()).find("TChiHH") != string::npos) ? true : false;
      // Check met out of 20%
      if (isSignal) {if (get<4>(it) < b.met()*0.8 || get<4>(it) > b.met()*1.2) continue;}
      else if (get<4>(it) != round(b.met())) continue;
      pair<int, int> signal_mass = getSignalMassValues(*b.FileNames().begin());
      // Check nlsp_mass
      if (signal_mass.first != get<5>(it)) continue;
      // Check lsp_mass
      int regionLsp = get<6>(it)==1? 0:get<6>(it);
      if (signal_mass.second != regionLsp) continue;
      // Check sample name
      //if (get<5>(it) != -1) cout<<*b.FileNames().begin()<<" "<<get<0>(it)<<endl;
      if ((*b.FileNames().begin()).find("TChiHH") != string::npos && get<0>(it).find("TChiHH") != string::npos)  {
        if ((*b.FileNames().begin()).find("HToBB_2D") != get<0>(it).find("2D")) continue;
      } else if (!findStringIC(*b.FileNames().begin(), get<0>(it))) continue;

      //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
      //cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<get<4>(it)<<endl;

      inRegion = true;
      break;
    }
  }
  return inRegion;
}

const NamedFunc boostSignalRegion("boostSignalRegion",[](const Baby &b) -> NamedFunc::ScalarType{
  return regionCut(b, 0);
});
const NamedFunc boostControlRegion("boostControlRegion",[](const Baby &b) -> NamedFunc::ScalarType{
  return regionCut(b, 1);
});

//const NamedFunc boostSignalRegion("boostSignalRegion",[](const Baby &b) -> NamedFunc::ScalarType{
//  bool inRegion = false;
//  int regionIndex = 0;
//  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int, int, int, int> > > > *> (b.EventVetoData());
//
//  // Print overlapping event numbers
//  for (auto const & itEventNumber : (*eventNumberData)[regionIndex]) {
//    if (itEventNumber.second.size() !=1) {
//      for (auto const & it : itEventNumber.second)
//        cout<<itEventNumber.first<<" "<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<get<4>(it)<<" "<<get<5>(it)<<endl;
//    }
//  }
//
//  if ((*eventNumberData)[regionIndex].count(b.event())==1) {
//    for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
//      //if (get<1>(it) == abs(b.SampleType()) && get<3>(it) == b.lumiblock() && get<2>(it) == b.run()) {
//      if (get<1>(it) == abs(b.SampleType()) && get<2>(it) == b.run()) {
//        if (findStringIC(*b.FileNames().begin(), get<0>(it))) {
//          //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
//          inRegion = true;
//          break;
//        }
//      }
//    }
//  }
//  return inRegion;
//});
//
//const NamedFunc boostControlRegion("boostControlRegion",[](const Baby &b) -> NamedFunc::ScalarType{
//  bool inRegion = false;
//  int regionIndex = 1;
//  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = static_cast<vector<map<Long64_t, set< tuple<string, int, int, int, int, int, int> > > > *> (b.EventVetoData());
//
//  // Print overlapping event numbers
//  for (auto const & itEventNumber : (*eventNumberData)[regionIndex]) {
//    if (itEventNumber.second.size() !=1) {
//      for (auto const & it : itEventNumber.second)
//        cout<<itEventNumber.first<<" "<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<get<4>(it)<<" "<<get<5>(it)<<endl;
//    }
//  }
//
//  if ((*eventNumberData)[regionIndex].count(b.event())==1) {
//    for (auto const & it : (*eventNumberData)[regionIndex][b.event()]) {
//      //if (get<1>(it) == abs(b.SampleType()) && get<3>(it) == b.lumiblock() && get<2>(it) == b.run()) {
//      if (get<1>(it) == abs(b.SampleType()) && get<2>(it) == b.run()) {
//        if (findStringIC(*b.FileNames().begin(), get<0>(it))) {
//          //if (b.met()<300) cout<<get<0>(it)<<" "<<get<1>(it)<<" "<<get<2>(it)<<" "<<get<3>(it)<<" "<<b.event()<<endl;
//          inRegion = true;
//          break;
//        }
//      }
//    }
//  }
//  return inRegion;
//});

string getLuminosityString(string const & year_string) {
  set<int> years;
  HigUtilities::parseYears(year_string, years);
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
  string const & nanoAodFolder, string const & production, string const & dataType, string const & fileTag, string const & sample_name, string const & year_string, 
  vector<shared_ptr<Process> > & procs) {
  set<int> years;
  //years_a = {2016};
  HigUtilities::parseYears(year_string, years);

  // Set base tags
  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
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
    if (fileTag=="all") fileNames = {"*TChiHH*.root"};
    else fileNames = {"*TChiHH_mChi-"+massPoints[0]+"_mLSP-"+massPoints[1]+"_*.root"};
  }
  else fileNames = mctags[fileTag];

  // Set folders
  torch::OrderedDict<string, string> folderDict;
  string mcProductionFolder = nanoAodFolder+"/"+production;
  string dataProductionFolder = nanoAodFolder+"/"+production;
  string signalProductionFolder = nanoAodFolder+"/"+production;
  folderDict.insert("mc_production_folder", mcProductionFolder);
  // preselect:
  //   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
  //   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
  //folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_higloose/");
  // higloose: 
  //   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
  //   met>150 && nvlep==0
  folderDict.insert("search_mc_skim_folder", "mc/merged_higmc_higloose/");
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
  folderDict.insert("ttbar_data_skim_folder","data/merged_higdata_higlep1T/");
  folderDict.insert("zll_data_skim_folder","data/merged_higdata_higlep2T/");
  folderDict.insert("qcd_data_skim_folder","data/merged_higdata_higqcd/");

  folderDict.insert("signal_production_folder", signalProductionFolder);
  folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higloose/");
  folderDict.insert("ttbar_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higlep1T/");
  folderDict.insert("zll_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higlep2T/");
  folderDict.insert("qcd_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higqcd/");

  // Set paths
  set<string> pathNames;
  pathNames = attach_folder(folderDict[dataType+"_production_folder"], years, folderDict[sample_name+"_"+dataType+"_skim_folder"], fileNames);

  if (dataType == "data") procs.push_back(Process::MakeShared<Baby_pico>(processName, type, color, pathNames, additionalCut));
  else procs.push_back(Process::MakeShared<Baby_pico>(processName, type, color, pathNames,"stitch"&&additionalCut)); // mc, signal
}

void addAllMcProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & sample_name, string const & year_string, 
  vector<shared_ptr<Process> > & procs) {
  Palette colors("txt/colors.txt", "default");
  addProcess("t#bar{t}+X (#tau_{had}>0) "+processName_postfix, Process::Type::background, colors("tt_htau"), "ntrutauh>0"&&additionalCut,
    nanoAodFolder, production, "mc", "tt", sample_name, year_string,
    procs);
  addProcess("t#bar{t}+X (#tau_{had}=0) "+processName_postfix, Process::Type::background, colors("tt_1l"), "ntrutauh==0"&&additionalCut,
    nanoAodFolder, production, "mc", "tt", sample_name, year_string,
    procs);
  addProcess("Z+Jets "+processName_postfix, Process::Type::background, kOrange+1, additionalCut,
    nanoAodFolder, production, "mc", "zjets", sample_name, year_string,
    procs);
  addProcess("W+Jets "+processName_postfix, Process::Type::background, kGreen+1, additionalCut,
    nanoAodFolder, production, "mc", "wjets", sample_name, year_string,
    procs);
  addProcess("Single t "+processName_postfix, Process::Type::background, colors("single_t"), additionalCut,
    nanoAodFolder, production, "mc", "single_t", sample_name, year_string,
    procs);
  addProcess("QCD "+processName_postfix, Process::Type::background, colors("other"), additionalCut,
    nanoAodFolder, production, "mc", "qcd", sample_name, year_string,
    procs);
  addProcess("Other "+processName_postfix, Process::Type::background, kGray+2, additionalCut,
    nanoAodFolder, production, "mc", "other", sample_name, year_string,
    procs);
}

void addMultipleSignalProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & sample_name, string const & year_string, 
  vector<shared_ptr<Process> > & procs) 
{
  vector<string> sigm = {"175", "500", "950"};
  for (unsigned isig(0); isig<sigm.size(); isig++){
    addProcess("TChiHH("+sigm[isig]+",0) "+processName_postfix, Process::Type::signal, 1, additionalCut,
      nanoAodFolder, production, "signal", sigm[isig]+"_0", sample_name, year_string,
      procs);
  }
}

void addDataProcess(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & sample, string const & year_string, 
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
      nanoAodFolder, production, "data", "", sample, year_string,
      procs);
}

//// Version with no MET
//void addEventNumberData(set<string> eventNumberFilenames, vector<map<Long64_t, set<tuple<string, int, int, int> > > > * eventNumberData) {
//  eventNumberData->push_back(map<Long64_t, set<tuple<string, int, int, int> > > ());
//  for (string const & eventNumberFilename : eventNumberFilenames) {
//    ifstream eventNumberFile(eventNumberFilename);
//    string line;
//    while (getline(eventNumberFile, line)) {
//      if(line[0]=='#') continue;
//      vector<string> lineSplit;
//      HigUtilities::stringToVectorString(line, lineSplit, ",");
//      string const & sampleType = lineSplit[0];
//      int year = stoi(lineSplit[1]);
//      int run = stoi(lineSplit[2]);
//      int lumiblock = stoi(lineSplit[3]);
//      Long64_t event = stoll(lineSplit[4]);
//      // Insert data
//      if ((*eventNumberData).back().count(event) == 0) {
//        (*eventNumberData).back()[event] = {{sampleType, year, run, lumiblock}};
//      } else {
//        if ((*eventNumberData).back()[event].count({sampleType, year, run, lumiblock}) != 0) cout<<"Duplicate: "<<sampleType<<" "<<year<<" "<<run<<" "<<lumiblock<<" "<<event<<endl;
//        (*eventNumberData).back()[event].insert({sampleType, year, run, lumiblock});
//      }
//    }
//    eventNumberFile.close();
//  }
//}

void addEventNumberData(set<string> eventNumberFilenames, vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData) {
  eventNumberData->push_back(map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > ());
  for (string const & eventNumberFilename : eventNumberFilenames) {
    ifstream eventNumberFile(eventNumberFilename);
    string line;
    while (getline(eventNumberFile, line)) {
      if(line[0]=='#') continue;
      if(line.find("tree") != string::npos) continue;
      vector<string> lineSplit;
      HigUtilities::stringToVectorString(line, lineSplit, ",");
      string const & sampleType = lineSplit[0];
      int year = stoi(lineSplit[1]);
      int run = stoi(lineSplit[2]);
      int lumiblock = stoi(lineSplit[3]);
      Long64_t event = stoll(lineSplit[4]);
      int met = stoi(lineSplit[5]);
      int nlsp_mass = -1; int lsp_mass = -1;
      bool isSignal = false; if (sampleType.find("TChiHH") != string::npos) isSignal = true;
      if (isSignal) {nlsp_mass = stoi(lineSplit[6]); lsp_mass = stoi(lineSplit[7]);}
      // Insert data
      if ((*eventNumberData).back().count(event) == 0) {
        (*eventNumberData).back()[event] = {{sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass}};
      } else {
        if ((*eventNumberData).back()[event].count({sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass}) != 0) cout<<"Duplicate: "<<sampleType<<" "<<year<<" "<<run<<" "<<lumiblock<<" "<<event<<" "<<met<<" "<<nlsp_mass<<" "<<lsp_mass<<endl;
        (*eventNumberData).back()[event].insert({sampleType, year, run, lumiblock, met, nlsp_mass, lsp_mass});
      }
    }
    eventNumberFile.close();
  }
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
  //string production_a = "higgsino_eldorado"; string nanoAodFolder_a = "/net/cms29/cms29r0/pico/NanoAODv5";
  //string higgsino_version = "v3";
  string higgsino_version = "";

  //string production_a = "higgsino_humboldt"+higgsino_version; string nanoAodFolder_a = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv5";
  string production_a = "higgsino_klamath"; 
  string nanoAodFolder_a = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7";
  string sample_a = sample_name;
  string year_string_a = year_string;

  string total_luminosity_string = HigUtilities::getLuminosityString(year_string_a);

  //    Define processes, including intersections
  //--------------------------------------------------
  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters =  Functions::hem_veto && HigUtilities::pass_2016;
  //NamedFunc base_filters = Functions::hem_veto && "pass && weight < 10";
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
  addAllMcProcesses("", "1",
    nanoAodFolder_a, production_a, 
    sample_a, year_string_a,
    procs);
  addMultipleSignalProcesses("", "1",
    nanoAodFolder_a, production_a, 
    sample_a, year_string_a,
    procs);
  if (unblind) {
    addDataProcess("", "1",
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs);
  }

  vector<shared_ptr<Process> > procs_all;
  addProcess("MC", Process::Type::background, 1, "1",
    nanoAodFolder_a, production_a, "mc", "all", sample_a, year_string_a,
    procs_all);
  addProcess("TChiHH", Process::Type::signal, 1, "1",
    nanoAodFolder_a, production_a, "signal", "all", sample_a, year_string_a,
    procs_all);
  if (unblind) {
    addDataProcess("DATA", "1",
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs_all);
  }

  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*Higfuncs::w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*Higfuncs::w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*Higfuncs::w_years*Functions::w_pileup;
  NamedFunc weight = Higfuncs::final_weight;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2_v0*Higfuncs::w_years;

  // Load event number data
  //set<string> vetoEventListNames = {"processed_resolved_list_CR_SCAN_MC.txt", "processed_resolved_list_SR_SCAN_MC.txt",
  //                                  "processed_resolved_list_CR_SCAN_TChiHH.txt", "processed_resolved_list_SR_SCAN_TChiHH.txt"};
  //set<string> vetoEventListNames = {"processed_resolved_list_CR_SCAN_MC.txt", "processed_resolved_list_SR_SCAN_MC.txt"};
  //set<string> boostedSignalRegion = {"BoostedEvents/boostedEvts_MC2016bkg_SR_noVeto.txt",
  //                                   "BoostedEvents/boostedEvts_MC2017bkg_SR_noVeto.txt",
  //                                   "BoostedEvents/boostedEvts_MC2018bkg_SR_noVeto.txt",
  //                                   "eventlist/processed_resolved_list_SR_SCAN_MC.txt"
  //                                   };
  //set<string> boostedControlRegion = {"BoostedEvents/boostedEvts_MC2016bkg_CR_noVeto.txt",
  //                                    "BoostedEvents/boostedEvts_MC2017bkg_CR_noVeto.txt",
  //                                    "BoostedEvents/boostedEvts_MC2018bkg_CR_noVeto.txt",
  //                                    "eventlist/processed_resolved_list_CR_SCAN_MC.txt"
  //                                    };
  set<string> boostedSignalRegion = {
                                     "BoostedEvents/boostedEvts_MC2016bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018bkg_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig_SR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig_SR_noVeto.txt",
                                     };
  set<string> boostedControlRegion = {
                                     "BoostedEvents/boostedEvts_MC2016bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018bkg_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2016sig_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2017sig_CR_noVeto.txt",
                                     "BoostedEvents/boostedEvts_MC2018sig_CR_noVeto.txt",
                                     };
  set<string> resolvedSignalRegion = {
                                      "eventlist/processed_resolved_list_SR_SCAN_MC.txt",
                                      "eventlist/processed_resolved_list_SR_SCAN_TChiHH1D.txt",
                                      //"eventlist/processed_resolved_list_SR_SCAN_TChiHH2D.txt",
                                     };
  set<string> resolvedControlRegion = {
                                      "eventlist/processed_resolved_list_CR_SCAN_MC.txt",
                                      "eventlist/processed_resolved_list_CR_SCAN_TChiHH1D.txt",
                                      //"eventlist/processed_resolved_list_CR_SCAN_TChiHH2D.txt",
                                      };
  //eventNumberData[0][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Boosted SR
  //eventNumberData[1][event] = {[sampleType, year, run, lumiblock, met, NSLP, LSP]} // Boosted CR
  //eventNumberData[2][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Resolved SR
  //eventNumberData[3][event] = {[sampleType, year, run, lumiblock, met, NLSP, LSP]} // Resolved CR
  vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > > * eventNumberData = new vector<map<Long64_t, set<tuple<string, int, int, int, int, int, int> > > >();
  addEventNumberData(boostedSignalRegion, eventNumberData);
  addEventNumberData(boostedControlRegion, eventNumberData);
  addEventNumberData(resolvedSignalRegion, eventNumberData);
  addEventNumberData(resolvedControlRegion, eventNumberData);

  //// Compare overlap between regions
  //{
  //  int regionA = 1;
  //  int regionB = 3;
  //  // (year, event), (sampleName, met, nlsp, lsp)
  //  // 0: sampleType, 1: year, 2: run, 3: lumiblock, 4: met, 5: nlsp_mass, 6: lsp_mass
  //  map<tuple<int, Long64_t>, set<tuple<string, int, int, int> > > eventDataA;
  //  cout<<"Fill eventData"<<endl;
  //  // Fill eventData for regionA
  //  for (auto const & itEventNumber : (*eventNumberData)[regionA]) {
  //    Long64_t const & event = itEventNumber.first;
  //    bool foundDuplicate = false;
  //    for (auto const & it : itEventNumber.second) {
  //      string const & sampleName = get<0>(it);
  //      int const & met = get<4>(it);
  //      int const & year = get<1>(it);
  //      int const nlsp_mass = get<5>(it);
  //      int const lsp_mass = get<6>(it);
  //      tuple<int, Long64_t> event_identifier = {year, event};
  //      tuple<string, int, int, int> event_data = {sampleName, met, nlsp_mass, lsp_mass};
  //      if (eventDataA.count(event_identifier) == 0) {
  //        eventDataA[event_identifier];
  //        eventDataA[event_identifier].insert(event_data);
  //      } else {
  //        if (eventDataA[event_identifier].count(event_data) == 0) eventDataA[event_identifier].insert(event_data);
  //        else {
  //          foundDuplicate = true;
  //          cout<<"[Duplicate] "<<sampleName<<" "<<year<<" "<<event<<" "<<met<<endl;
  //        }
  //      }
  //      //if (met<300) foundDuplicate = true;
  //    }
  //    if (foundDuplicate) {
  //      for (auto const & itPrint : itEventNumber.second) {
  //        cout<<regionA<<" "<<itEventNumber.first<<" "<<get<0>(itPrint)<<" "<<get<1>(itPrint)<<" "<<get<2>(itPrint)<<" "<<get<3>(itPrint)<<" "<<get<4>(itPrint)<<" "<<get<5>(itPrint)<<" "<<get<6>(itPrint)<<endl;
  //      }
  //    }
  //  }
  //  // (year, event), (sampleName, met, nlsp, lsp)
  //  map<tuple<int, Long64_t>, set<tuple<string, int, int, int> > > eventDataB;
  //  cout<<"Fill eventData"<<endl;
  //  // Fill eventData for regionB
  //  for (auto const & itEventNumber : (*eventNumberData)[regionB]) {
  //    Long64_t const & event = itEventNumber.first;
  //    bool foundDuplicate = false;
  //    for (auto const & it : itEventNumber.second) {
  //      string const & sampleName = get<0>(it);
  //      int const & met = get<4>(it);
  //      int const & year = get<1>(it);
  //      int const nlsp_mass = get<5>(it);
  //      int const lsp_mass = get<6>(it);
  //      tuple<int, Long64_t> event_identifier = {year, event};
  //      tuple<string, int, int, int> event_data = {sampleName, met, nlsp_mass, lsp_mass};
  //      if (eventDataB.count(event_identifier) == 0) {
  //        eventDataB[event_identifier];
  //        eventDataB[event_identifier].insert(event_data);
  //      } else {
  //        if (eventDataB[event_identifier].count(event_data) == 0) eventDataB[event_identifier].insert(event_data);
  //        else {
  //          foundDuplicate = true;
  //          cout<<"[Duplicate] "<<sampleName<<" "<<year<<" "<<event<<" "<<met<<endl;
  //        }
  //      }
  //    }
  //    if (foundDuplicate) {
  //      for (auto const & itPrint : itEventNumber.second) {
  //        cout<<regionB<<" "<<itEventNumber.first<<" "<<get<0>(itPrint)<<" "<<get<1>(itPrint)<<" "<<get<2>(itPrint)<<" "<<get<3>(itPrint)<<" "<<get<4>(itPrint)<<" "<<get<5>(itPrint)<<" "<<get<6>(itPrint)<<endl;
  //      }
  //    }
  //  }

  //  // Find overlap between regionA and regionB
  //  for (auto const & it : eventDataA) {
  //    if (eventDataB.count(it.first) != 0) {
  //      auto const & set_eventA = eventDataA[it.first];
  //      auto const & set_eventB = eventDataB[it.first];
  //      // Check nlsp and lsp
  //      bool sameNlsp = false;
  //      for (auto const & eventDetailA: set_eventA) {
  //        for (auto const & eventDetailB: set_eventB) {
  //          if (get<2>(eventDetailA) == get<2>(eventDetailB)) {
  //            //// Find events with met larger than 80% different
  //            //if (get<1>(eventDetailA) < get<1>(eventDetailB)*0.8 || get<1>(eventDetailA) > get<1>(eventDetailB)*1.2) sameNlsp = true;
  //            if (get<1>(eventDetailA) > get<1>(eventDetailB)*0.7 && get<1>(eventDetailA) < get<1>(eventDetailB)*1.3) sameNlsp = true;
  //          }
  //        }
  //      }
  //      if (sameNlsp) {
  //        cout<<"[Common event] "<<get<0>(it.first)<<" "<<get<1>(it.first)<<endl;
  //        for (auto const & itSet : set_eventA) cout<<regionA<<" "<<get<0>(itSet)<<" "<<get<1>(itSet)<<" "<<get<2>(itSet)<<" "<<get<3>(itSet)<<endl;
  //        for (auto const & itSet : set_eventB) cout<<regionB<<" "<<get<0>(itSet)<<" "<<get<1>(itSet)<<" "<<get<2>(itSet)<<" "<<get<3>(itSet)<<endl;
  //      }
  //    }
  //  }
  //}

  //// 0: sampleName, 1: year, 2: run, 3: lumiblock, 4: met, 5: nlsp_mass, 6: lsp_mass
  //// Print overlapping event numbers with same (sampleName, year, met)
  //for (unsigned regionIndex = 0; regionIndex < eventNumberData->size(); ++regionIndex) {
  //  for (auto const & itEventNumber : (*eventNumberData)[regionIndex]) {
  //    if (itEventNumber.second.size() !=1) {
  //      set<tuple<string, int, int> > sample;
  //      bool foundDuplicate = false;
  //      for (auto const & it : itEventNumber.second) {
  //        string const & sampleName = get<0>(it);
  //        int const & met = get<4>(it);
  //        int const & year = get<1>(it);
  //        if (sample.count({sampleName, year, met}) == 0) {
  //          sample.insert({sampleName, year, met});
  //        } else {
  //          for (auto const & itSample : sample) cout<<get<0>(itSample)<<" "<<get<1>(itSample)<<" "<<get<2>(itSample)<<endl;
  //          foundDuplicate = true;
  //          cout<<itEventNumber.first<<" "<<sampleName<<" "<<year<<" "<<met<<endl;
  //          break;
  //        }
  //      }
  //      if (foundDuplicate) {
  //        for (auto const & itPrint : itEventNumber.second) {
  //          cout<<regionIndex<<" "<<itEventNumber.first<<" "<<get<0>(itPrint)<<" "<<get<1>(itPrint)<<" "<<get<2>(itPrint)<<" "<<get<3>(itPrint)<<" "<<get<4>(itPrint)<<" "<<get<5>(itPrint)<<" "<<get<6>(itPrint)<<endl;
  //        }
  //      }
  //    }
  //  }
  //}

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


  NamedFunc signalRegion = "hig_cand_am[0]>100 && hig_cand_am[0]<140 && nbm>=3";

  NamedFunc rsr = base_filters&&resolved_cuts&&signalRegion;
  NamedFunc rcr = base_filters&&resolved_cuts&&!signalRegion;
  NamedFunc bsr = boostSignalRegion;
  NamedFunc bcr = boostControlRegion;

  //        Define and fill tables
  //------------------------------------
  PlotMaker pm;
  vector<TableRow> tablerows;
  for (auto &imet: vc_met) {
    // Label for bin
    tablerows.push_back(TableRow("$"+CodeToLatex((imet).Data())+"$"));

    NamedFunc rsr_bin = rsr&&imet;
    NamedFunc rcr_bin = rcr&&imet;
    NamedFunc bsr_bin = imet&&bsr;
    NamedFunc bcr_bin = imet&&bcr;

    // Add rows to table
    tablerows.push_back(TableRow("RSR", rsr_bin, 0,0, weight));
    tablerows.push_back(TableRow("RCR", rcr_bin, 0,0, weight));
    tablerows.push_back(TableRow("BSR", bsr_bin, 0,0, weight));
    tablerows.push_back(TableRow("BCR", bcr_bin, 0,0, weight));

    tablerows.push_back(TableRow("Exclusive Binning"));
    tablerows.push_back(TableRow("Pure RSR", rsr_bin && !(bcr_bin||bsr_bin), 0,0, weight));
    tablerows.push_back(TableRow("Pure RCR", rcr_bin && !(bcr_bin||bsr_bin), 0,0, weight));
    tablerows.push_back(TableRow("Pure BSR", bsr_bin && !(rcr_bin||rsr_bin), 0,0, weight));
    tablerows.push_back(TableRow("Pure BCR", bcr_bin && !(rcr_bin||rsr_bin), 0,1, weight));
    tablerows.push_back(TableRow("RSR and BSR", rsr_bin && bsr_bin, 0,0, weight));
    tablerows.push_back(TableRow("RSR and BCR", rsr_bin && bcr_bin, 0,0, weight));
    tablerows.push_back(TableRow("RCR and BSR", rcr_bin && bsr_bin, 0,0, weight));
    tablerows.push_back(TableRow("RCR and BCR", rcr_bin && bcr_bin, 0,0, weight));

  }
  TString tabname = "table_"; tabname += "resolved";

  bool print_uncertainty = true;
  pm.Push<Table>("FixName:table_"+sample_a+"_regions", tablerows, procs, 0, 1, 0, 1, 0, print_uncertainty).LuminosityTag(total_luminosity_string);

  //pm.Push<EventScan>("resolved_list_SR", base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event"}, procs_all, 10, true);
  //pm.Push<EventScan>("resolved_list_CR", base_filters&&search_resolved_cuts&&!signalRegion,vector<NamedFunc>{"SampleType", "run","lumiblock","event"}, procs_all, 10, true);
  //pm.Push<EventScan>("resolved_list_SR", base_filters&&search_resolved_cuts&&signalRegion&&eventNumberVeto,vector<NamedFunc>{"SampleType","run","lumiblock","event"}, procs_all, 10, true);
  //pm.Push<EventScan>("resolved_list_CR", base_filters&&search_resolved_cuts&&!signalRegion&&eventNumberVeto,vector<NamedFunc>{"SampleType", "run","lumiblock","event"}, procs_all, 10, true);
  //pm.Push<EventScan>("resolved_list_SR_met300", "met>300"&&base_filters&&search_resolved_cuts&&signalRegion,vector<NamedFunc>{"SampleType","run","lumiblock","event", "met", "nbl", "nbm", "nbt", "hig_cand_am[0]","hig_cand_dm[0]", "hig_cand_drmax[0]"}, procs_all, 10, true);
  //pm.Push<EventScan>("resolved_list_CR_met300", "met>300"&&base_filters&&search_resolved_cuts&&!signalRegion,vector<NamedFunc>{"SampleType", "run","lumiblock","event", "met", "nbl", "nbm", "nbt", "hig_cand_am[0]","hig_cand_dm[0]", "hig_cand_drmax[0]"}, procs_all, 10, true);

  pm.SetEventVetoData(static_cast<void * >(eventNumberData));

  pm.min_print_ = true;
  pm.MakePlots(1.);

  //// Make new event list with processed_ prefix.
  ////set<string> eventListNames = {"resolved_list_CR_SCAN_MC.txt","resolved_list_SR_SCAN_MC.txt",
  ////                              "resolved_list_CR_SCAN_TChiHH.txt", "resolved_list_SR_SCAN_TChiHH.txt"};
  //set<string> eventListNames = {"resolved_list_CR_met300_SCAN_MC.txt","resolved_list_SR_met300_SCAN_MC.txt",
  //                              "resolved_list_CR_met300_SCAN_TChiHH.txt", "resolved_list_SR_met300_SCAN_TChiHH.txt"};
  //for (string const & eventListName : eventListNames) {
  //  ifstream eventList(eventListName);
  //  string processedEventListName = "processed_"+eventListName;
  //  cout<<"Processing "<<eventListName<<" to "<<processedEventListName<<endl;
  //  ofstream processedEventList(processedEventListName);
  //  processedEventList<<"# SampleName, Year, RunNumber, LumiBlockNumber, EventNumber"<<endl;
  //  string line;
  //  while (getline(eventList, line)) {
  //    if (Contains(line, "Instance")) continue;
  //    vector<string> lineSplit;
  //    HigUtilities::stringToVectorString(line, lineSplit, " ");
  //    // Expect "Row Instance SampleType run lumiblock event Filename" sequence
  //    int sampleType = stoi(lineSplit[2]); string run = lineSplit[3]; string lumiblock = lineSplit[4]; string event = lineSplit[5]; string filename = lineSplit[6];
  //    string datasetName = getPartialDatasetName(filename);
  //    processedEventList<<datasetName<<", "<<abs(sampleType)<<", "<<run<<", "<<lumiblock<<", "<<event<<endl;
  //  }
  //  eventList.close();
  //  processedEventList.close();
  //};

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 0},
      {"year", required_argument, 0, 0},
      {"unblind", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "nl:t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 0:
      optname = long_options[option_index].name;
      if(optname == "sample"){
        sample_name = optarg;
      } else if (optname == "year") {
        year_string = optarg;
      } else if (optname == "unblind") {
        unblind = true;
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
