///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <getopt.h>

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

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  //if (b.type()==106000) {
  //  return 35.9;
  //}
  if (b.SampleType()==2016){
    return weight*35.9;
  } else if (b.SampleType()==2017){
    return weight*41.5;
  } else {
    return weight*59.6;
  }
});

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
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  set<string> fileNames;
  if (dataType == "data") fileNames = {"*.root"};
  else if (dataType == "signal") {
    vector<string> massPoints;
    HigUtilities::stringToVectorString(fileTag, massPoints, "_");
    // Special case for angeles production
    if (production == "higgsino_angeles" && massPoints[1] == "0") massPoints[1] = "1";
    fileNames = {"*TChiHH_mChi-"+massPoints[0]+"_mLSP-"+massPoints[1]+"_*.root"};
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
  //   (Max$(el_pt*el_sig)>40 || Max$(mu_pt*mu_sig)>40) // pass_1l_trig40 (sig is signal lepton)
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
    NamedFunc triggers_data = "1";
    NamedFunc lepton_triggers = "(HLT_IsoMu24 || HLT_IsoMu27 || HLT_Mu50 || HLT_Ele27_WPTight_Gsf || HLT_Ele35_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT)";
    NamedFunc met_triggers = "(HLT_PFMET110_PFMHT110_IDTight || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight || HLT_PFMET120_PFMHT120_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFMET120_PFMHT120_IDTight_PFHT60 || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60)";
    if (sample == "zll") triggers_data = lepton_triggers;
    else if (sample == "ttbar") triggers_data = lepton_triggers || met_triggers;
    else if (sample == "qcd") triggers_data = met_triggers;
    else  triggers_data = met_triggers;

    addProcess("Data"+processName_postfix, Process::Type::data, kBlack, triggers_data && additionalCut,
      nanoAodFolder, production, "data", "", sample, year_string,
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
  //string production_a = "higgsino_eldorado"; string nanoAodFolder_a = "/net/cms29/cms29r0/pico/NanoAODv5";
  string production_a = "higgsino_humboldt"; string nanoAodFolder_a = "/net/cms25/cms25r5/pico/NanoAODv5";
  string sample_a = sample_name;
  string year_string_a = year_string;

  string total_luminosity_string = getLuminosityString(year_string_a);

  //    Define processes, including intersections
  //--------------------------------------------------
  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...
  NamedFunc base_filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  NamedFunc search_resolved_cuts = 
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved_cuts = 
                         "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc zll_resolved_cuts =
                         "nlep==2&&njet>=4&&njet<=5&&met<50&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  NamedFunc qcd_resolved_cuts =
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";

  NamedFunc resolved_cuts = "1";
  if (sample_a == "ttbar") resolved_cuts = ttbar_resolved_cuts;
  else if (sample_a == "zll") resolved_cuts = zll_resolved_cuts;
  else if (sample_a == "qcd") resolved_cuts = qcd_resolved_cuts;
  else resolved_cuts = search_resolved_cuts;

  vector<shared_ptr<Process> > procs;
  addAllMcProcesses("", base_filters&&resolved_cuts,
    nanoAodFolder_a, production_a, 
    sample_a, year_string_a,
    procs);
  addMultipleSignalProcesses("", base_filters&&resolved_cuts,
    nanoAodFolder_a, production_a, 
    sample_a, year_string_a,
    procs);
  if (unblind) {
    addDataProcess("", base_filters&&resolved_cuts,
      nanoAodFolder_a, production_a, 
      sample_a, year_string_a,
      procs);
  }

  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig_run2*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig*w_years;
  NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;


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

  pm.min_print_ = true;
  pm.MakePlots(1.);

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
