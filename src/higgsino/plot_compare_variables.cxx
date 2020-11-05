#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/ordered_dict.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
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

// nanoAodFolder: /net/cms29/cms29r0/pico/NanoAODv5/
// production: higgsino_humboldt
// dataType: mc/data/signal
// fileTag: (tt, single_t, zjets, wjets, qcd, other, all), (all), (200_1, 600_1, 950_1, ...)
// sample: search/ttbar/zll/qcd
// year_string: 2016/2017/2018/run2
void addProcess(string const & processName, Process::Type type, int color, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & dataType, string const & fileTag, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs) {
  set<int> years;
  //years_a = {2016};
  if (year_string == "run2") years = {2016, 2017, 2018};
  else HigUtilities::parseYears(year_string, years);

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
  folderDict.insert("search_data_skim_folder","data/merged_higmc_higloose/");
  folderDict.insert("ttbar_data_skim_folder","data/merged_higmc_higlep1T/");
  folderDict.insert("zll_data_skim_folder","data/merged_higmc_higlep2T/");
  folderDict.insert("qcd_data_skim_folder","data/merged_higmc_higqcd/");

  folderDict.insert("signal_production_folder", signalProductionFolder);
  folderDict.insert("search_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higloose/");
  folderDict.insert("ttbar_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higlep1T/");
  folderDict.insert("zll_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higlep2T/");
  folderDict.insert("qcd_signal_skim_folder", "SMS-TChiHH_2D/merged_higmc_higqcd/");

  if (sample == "ttbar") folderDict.insert("mc_skim_folder", folderDict["ttbar_mc_skim_folder"]);
  else if (sample == "zll") folderDict.insert("mc_skim_folder", folderDict["zll_mc_skim_folder"]);
  else if (sample == "qcd") folderDict.insert("mc_skim_folder", folderDict["qcd_mc_skim_folder"]);
  else folderDict.insert("mc_skim_folder", folderDict["search_mc_skim_folder"]);

  // Set paths
  set<string> pathNames;
  pathNames = attach_folder(folderDict[dataType+"_production_folder"], years, folderDict[sample+"_"+dataType+"_skim_folder"], fileNames);

  procs.push_back(Process::MakeShared<Baby_pico>(processName, type, color, pathNames,"stitch"&&additionalCut));
}

void addProcess(string const & processName, int color, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & dataType, string const & fileTag, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs) 
{
  // set Type
  Process::Type type;
  if (dataType == "data") {
    type = Process::Type::data;
  } else if (dataType == "signal") {
    type = Process::Type::signal;
  } else { // mc
    type = Process::Type::background;
  }
  addProcess(processName, type, color, additionalCut,
    nanoAodFolder, production, dataType, fileTag, sample, year_string,
    procs);
}

void addProcess(string const & processName, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & dataType, string const & fileTag, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs) 
{
  Palette colors("txt/colors.txt", "default");

  // set Type and color
  Process::Type type;
  int color = 0;
  if (dataType == "data") {
    type = Process::Type::data;
    color = kBlack;
  } else if (dataType == "signal") {
    type = Process::Type::signal;
    // TODO determine color better with fileTag
    if (fileTag == "200_1") color = kGreen+1;
    else if (fileTag == "600_1") color = kRed;
    else color = kBlue;
  } else { // mc
    type = Process::Type::background;
    if (fileTag == "single_t") color = colors("single_t");
    else if (fileTag == "zjets") color = kOrange+1;
    else if (fileTag == "wjets") color = kGreen+1;
    else if (fileTag == "qcd") color = colors("other");
    else if (fileTag == "other") color = kGray+2;
    else { // tt
      if (Contains(additionalCut.Name(),"ntrutauh>0")) color = colors("tt_htau");
      else color = colors("tt_1l");
    }
  }
  addProcess(processName, type, color, additionalCut,
    nanoAodFolder, production, dataType, fileTag, sample, year_string,
    procs);
}

void addAllMcProcesses(string const & processName_postfix, NamedFunc const & additionalCut, 
  string const & nanoAodFolder, string const & production, string const & dataType, string const & sample, string const & year_string, 
  vector<shared_ptr<Process> > & procs) {
  Palette colors("txt/colors.txt", "default");
  addProcess("t#bar{t}+X (#tau_{had}>0) "+processName_postfix, Process::Type::background, colors("tt_htau"), "ntrutauh>0"&&additionalCut,
    nanoAodFolder, production, dataType, "tt", sample, year_string,
    procs);
  addProcess("t#bar{t}+X (#tau_{had}=0) "+processName_postfix, Process::Type::background, colors("tt_1l"), "ntrutauh==0"&&additionalCut,
    nanoAodFolder, production, dataType, "tt", sample, year_string,
    procs);
  addProcess("Z+Jets "+processName_postfix, Process::Type::background, kOrange+1, additionalCut,
    nanoAodFolder, production, dataType, "zjets", sample, year_string,
    procs);
  addProcess("W+Jets "+processName_postfix, Process::Type::background, kGreen+1, additionalCut,
    nanoAodFolder, production, dataType, "wjets", sample, year_string,
    procs);
  addProcess("Single t "+processName_postfix, Process::Type::background, colors("single_t"), additionalCut,
    nanoAodFolder, production, dataType, "single_t", sample, year_string,
    procs);
  addProcess("QCD "+processName_postfix, Process::Type::background, colors("other"), additionalCut,
    nanoAodFolder, production, dataType, "qcd", sample, year_string,
    procs);
  addProcess("Other "+processName_postfix, Process::Type::background, kGray+2, additionalCut,
    nanoAodFolder, production, dataType, "other", sample, year_string,
    procs);
}

void getSampleCutDict(string const & sample, torch::OrderedDict<string, NamedFunc> & sampleCutDict) {
  //NamedFunc base_resolved = 
  //                       "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  //NamedFunc ttbar_resolved_cuts = 
  //                       "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  //NamedFunc zll_resolved_cuts =
  //                       "nlep==2&&njet>=4&&njet<=5&&met<50&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  //NamedFunc qcd_resolved_cuts =
  //                       "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  torch::OrderedDict<string, NamedFunc> search_resolved_cuts;
  search_resolved_cuts.insert("ntk", "ntk==0");
  search_resolved_cuts.insert("low_dphi_met", "!low_dphi_met");
  search_resolved_cuts.insert("nvlep", "nvlep==0");
  search_resolved_cuts.insert("met", "met>150");
  search_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  search_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  search_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  search_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  search_resolved_cuts.insert("btags", "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))");
  search_resolved_cuts.insert("mht_filter", "met/mht<2");
  search_resolved_cuts.insert("met_calo_filter", "met/met_calo<2");
  search_resolved_cuts.insert("weight", "weight<10");
  torch::OrderedDict<string, NamedFunc> ttbar_resolved_cuts;
  ttbar_resolved_cuts.insert("nlep", "nlep==1");
  ttbar_resolved_cuts.insert("lep_pt", Higfuncs::lead_signal_lepton_pt>30);
  ttbar_resolved_cuts.insert("mt", "mt<=100");
  ttbar_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  ttbar_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  ttbar_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  ttbar_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  ttbar_resolved_cuts.insert("btags", "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))");
  ttbar_resolved_cuts.insert("met_calo_filter", "met/met_calo<5");
  ttbar_resolved_cuts.insert("weight", "weight<10");
  torch::OrderedDict<string, NamedFunc> zll_resolved_cuts;
  zll_resolved_cuts.insert("nlep", "nlep==2");
  zll_resolved_cuts.insert("lep_pt", Higfuncs::lead_signal_lepton_pt>30);
  zll_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  zll_resolved_cuts.insert("met", "met<50");
  zll_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  zll_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  zll_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  zll_resolved_cuts.insert("dbtags", "(nbm==0||nbm==1||nbm==2||nbm>=3)");
  zll_resolved_cuts.insert("met_calo_filter", "met/met_calo<5");
  zll_resolved_cuts.insert("weight", "weight<10");
  torch::OrderedDict<string, NamedFunc> qcd_resolved_cuts;
  qcd_resolved_cuts.insert("low_dphi_met", "low_dphi_met");
  qcd_resolved_cuts.insert("nvlep", "nvlep==0");
  qcd_resolved_cuts.insert("met", "met>150");
  qcd_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  qcd_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  qcd_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  qcd_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  qcd_resolved_cuts.insert("dbtags", "(nbm==0||nbm==1||nbm==2||nbm>=3)");
  qcd_resolved_cuts.insert("mht_filter", "met/mht<2");
  qcd_resolved_cuts.insert("met_calo_filter", "met/met_calo<2");
  qcd_resolved_cuts.insert("weight", "weight<10");


  if (sample == "ttbar") sampleCutDict = ttbar_resolved_cuts;
  else if (sample == "zll") sampleCutDict = zll_resolved_cuts;
  else if (sample == "qcd") sampleCutDict = qcd_resolved_cuts;
  else sampleCutDict = search_resolved_cuts; // search
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt log_shapes = log_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt log_shapes_info = log_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_log_shapes = {log_shapes};
  vector<PlotOpt> plt_log_shapes_info = {log_shapes_info};

  NamedFunc studyCut = "1";

  // Compare between productions or years or (search, control) samples or mc/data/signal
  NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years*Functions::w_pileup;
  //NamedFunc weight = "weight";
  //NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;

  // [Past option] Compares between productions, years, mc/data/signal.
  // [Current option] Compare between productions or years or (search, control) samples or mc/data/signal
  //   - Make different process for each n-1.
  //string sample = "search";

  // dataType: mc/data/signal
  // fileTag: (tt, single_t, zjets, wjets, qcd, other, all), (all), (200_1, 600_1, 950_1, ...)
  // sample: search/ttbar/zll/qcd
  // year_string: 2016/2017/2018/run2
  //string production_a = "higgsino_eldorado"; string nanoAodFolder_a = "/net/cms29/cms29r0/pico/NanoAODv5";
  string production_a = "higgsino_humboldt"; string nanoAodFolder_a = "/net/cms25/cms25r5/pico/NanoAODv5";
  string dataType_a = "mc";
  string fileTag_a = "all";
  string sample_a = "qcd";
  string year_string_a = "2016";

  //string production_b = "higgsino_eldorado"; string nanoAodFolder_b = "/net/cms29/cms29r0/pico/NanoAODv5";
  string production_b = "higgsino_humboldt"; string nanoAodFolder_b = "/net/cms25/cms25r5/pico/NanoAODv5";
  string dataType_b = "mc";
  string fileTag_b = "all";
  string sample_b = "qcd";
  string year_string_b = "2017";

  //NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2" && studyCut; //since pass_fsjets is not quite usable...
  NamedFunc base_filters = Functions::hem_veto && "pass";

  // Get sample cut dict: sampleCutDict["variable"] = "variable cut"
  torch::OrderedDict<string, NamedFunc> sampleCutDict_a;
  getSampleCutDict(sample_a, sampleCutDict_a);
  torch::OrderedDict<string, NamedFunc> sampleCutDict_b;
  getSampleCutDict(sample_b, sampleCutDict_b);

  torch::OrderedDict<string, Axis> axis_dict;
  axis_dict.insert("ntk",Axis(10, -0.5, 9.5, "ntk", "Number of iso tk", {0.5}));
  axis_dict.insert("npv",Axis(81, -0.5, 80.5, "npv", "Number of reconstructed primary vertices", {}));
  axis_dict.insert("npu_tru_mean",Axis(81, -0.5, 80.5, "npu_tru_mean", "Mean true interaction vertices", {}));
  axis_dict.insert("low_dphi_met",Axis(2, -0.5, 1.5, "low_dphi_met", "Low #Delta #phi", {0.5}));
  axis_dict.insert("nvlep",Axis(6, -0.5, 5.5, "nvlep", "N_{veto leps}", {0.5}));
  axis_dict.insert("nlep",Axis(5, 0.5, 5.5, "nlep", "N_{leps}", {}));
  axis_dict.insert("nmu",Axis(5, 0.5, 5.5, "nmu", "N_{muons}", {}));
  axis_dict.insert("nel",Axis(5, 0.5, 5.5, "nel", "N_{electrons}", {}));
  axis_dict.insert("met",Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}));
  axis_dict.insert("met_zll",Axis(20, 0, 150., "met", "p_{T}^{miss} [GeV]", {30}));
  axis_dict.insert("met_ttbar",Axis(16, 0, 800., "met", "p_{T}^{miss} [GeV]", {}));
  axis_dict.insert("ht",Axis(20, 0, 1000., "ht", "HT [GeV]", {}));
  axis_dict.insert("njet",Axis(12, -0.5, 11.5, "njet", "N_{jets}", {3.5, 5.5}));
  axis_dict.insert("hig_cand_drmax",Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}));
  axis_dict.insert("hig_cand_am",Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}));
  axis_dict.insert("hig_cand_dm",Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}));
  axis_dict.insert("btags",Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}));
  axis_dict.insert("lep_pt",Axis(10, 0, 300., "lep_pt[0]", "p_{l} [GeV]", {30}));
  axis_dict.insert("mt",Axis(10, 0, 200., "mt", "m_{T} [GeV]", {100}));
  axis_dict.insert("dbtags",Axis(5, -0.5, 4.5, "nbm", "N_{b medium}", {}));
  axis_dict.insert("ll_pt",Axis(16, 0, 400., "ll_pt[0]", "p_{T}^{ll} [GeV]", {75., 150., 200., 300}));
  axis_dict.insert("ll_m",Axis(16, 0, 400., "ll_m", "m^{ll} [GeV]", {80., 100.}));
  axis_dict.insert("mht_filter",Axis(10, 0, 10, "met/mht", "MET/MHT", {2}));
  axis_dict.insert("met_calo_filter",Axis(10, 0, 10, "met/met_calo", "MET/MET_calo", {2,5}));
  axis_dict.insert("weight",Axis(40, 0, 20, "weight", "weight", {10}));

  vector<string> log_plots = {"met", "btags"};

  // Make plots
  PlotMaker pm;

  set<string> drawVariables;
  // Find common variables to plot
  for (auto const & item : sampleCutDict_a) {
    if (sampleCutDict_b.find(item.key()) == nullptr) continue;
    drawVariables.insert(item.key());
  }
  // Add variables for special cases
  drawVariables.insert("ht");
  if (sample_a == "ttbar" || sample_b == "ttbar") drawVariables.insert("met");
  if (sample_a == "zll" || sample_b == "zll") drawVariables.insert("ll_pt");

  // Draw n-1
  // Loop over variables
  for (auto const & target_var : drawVariables) {
    // Make n-1 cut for target_var
    NamedFunc n_minus_1_cut_a = "1";
    for (auto const & cut_item : sampleCutDict_a) {
      if (cut_item.key() == target_var) continue;
      // Below are special conditions
      if (target_var == "njet") {
        if (cut_item.key() == "hig_cand_drmax") continue;
        if (cut_item.key() == "hig_cand_am") continue;
        if (cut_item.key() == "hig_cand_dm") continue;
      }
      if (target_var == "nlep") {
        if (cut_item.key() == "mt") continue;
      }
      n_minus_1_cut_a = n_minus_1_cut_a && cut_item.value();
    }
    // Make n-1 cut for target_var
    NamedFunc n_minus_1_cut_b = "1";
    for (auto const & cut_item : sampleCutDict_b) {
      if (cut_item.key() == target_var) continue;
      // Below are special conditions
      if (target_var == "njet") {
        if (cut_item.key() == "hig_cand_drmax") continue;
        if (cut_item.key() == "hig_cand_am") continue;
        if (cut_item.key() == "hig_cand_dm") continue;
      }
      if (target_var == "nlep") {
        if (cut_item.key() == "mt") continue;
      }
      n_minus_1_cut_b = n_minus_1_cut_b && cut_item.value();
    }

    vector<shared_ptr<Process> > procs;
    addProcess(CopyReplaceAll(production_a,"higgsino_","")+"_"+dataType_a+"_"+fileTag_a+"_"+sample_a+"_"+year_string_a, kGreen, n_minus_1_cut_a,
      nanoAodFolder_a, production_a, 
      dataType_a, fileTag_a, sample_a, year_string_a,
      procs);
    addProcess(CopyReplaceAll(production_b,"higgsino_","")+"_"+dataType_b+"_"+fileTag_b+"_"+sample_b+"_"+year_string_b, kBlue, n_minus_1_cut_b,
      nanoAodFolder_b, production_b, 
      dataType_b, fileTag_b, sample_b, year_string_b,
      procs);

    // Draw target_var where a and b are compared
    string figureName = "fig_n-1_"+target_var+"__"+production_a+"_"+dataType_a+"_"+fileTag_a+"_"+sample_a+"_"+CopyReplaceAll(year_string_a, ",","_")+"_vs_"+production_b+"_"+dataType_b+"_"+fileTag_b+"_"+sample_b+"_"+CopyReplaceAll(year_string_b, ",","_");
    // For special case
    if (target_var == "met" && (Contains(sample_a, "zll") || Contains(sample_b, "zll"))) {
      pm.Push<Hist1D>(axis_dict[target_var+"_zll"],
        base_filters,
        procs, plt_shapes).Weight(weight).Tag("FixName:"+figureName);
    } else if (target_var == "met" && (Contains(sample_a, "ttbar") || Contains(sample_b, "ttbar"))) {
      pm.Push<Hist1D>(axis_dict[target_var+"_ttbar"],
        base_filters,
        procs, plt_shapes).Weight(weight).Tag("FixName:"+figureName);
    // For log plots
    } else if (std::find(log_plots.begin(), log_plots.end(), target_var) != log_plots.end()) {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters,
        procs, plt_log_shapes).Weight(weight).Tag("FixName:"+figureName);
    // Normal case
    } else {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters,
        procs, plt_shapes).Weight(weight).Tag("FixName:"+figureName);
    }

    // In case dataType is mc with all, draw variables separately for mc
    if (dataType_a == "mc" && fileTag_a == "all") {
      vector<shared_ptr<Process> > procs_mc;
      addAllMcProcesses(CopyReplaceAll(production_a,"higgsino_","")+"_"+year_string_a, "1",
        nanoAodFolder_a, production_a, 
        dataType_a, sample_a, year_string_a,
        procs_mc);
      // Draw target_var
      string figureName_a = "fig_n-1_"+target_var+"__"+production_a+"_"+dataType_a+"_"+fileTag_a+"_"+sample_a+"_"+CopyReplaceAll(year_string_a, ",","_");
      if (target_var == "met" && Contains(sample_a, "zll")) {
        pm.Push<Hist1D>(axis_dict[target_var+"_zll"],
          base_filters && n_minus_1_cut_a,
          procs_mc, plt_lin).Weight(weight).Tag("FixName:"+figureName_a);
      } else if (target_var == "met" && Contains(sample_a, "ttbar")) {
        pm.Push<Hist1D>(axis_dict[target_var+"_ttbar"],
          base_filters && n_minus_1_cut_a,
          procs_mc, plt_log).Weight(weight).Tag("FixName:"+figureName_a);
      // For log plots
      } else if (std::find(log_plots.begin(), log_plots.end(), target_var) != log_plots.end()) {
        pm.Push<Hist1D>(axis_dict[target_var],
          base_filters && n_minus_1_cut_a,
          procs_mc, plt_log).Weight(weight).Tag("FixName:"+figureName_a);
      // Normal case
      } else {
        pm.Push<Hist1D>(axis_dict[target_var],
          base_filters && n_minus_1_cut_a,
          procs_mc, plt_lin).Weight(weight).Tag("FixName:"+figureName_a);
      }
    }

    if (dataType_b == "mc" && fileTag_b == "all") {
      vector<shared_ptr<Process> > procs_mc;
      addAllMcProcesses(CopyReplaceAll(production_b,"higgsino_","")+"_"+year_string_b, "1",
        nanoAodFolder_b, production_b, 
        dataType_b, sample_b, year_string_b,
        procs_mc);
      // Draw target_var
      string figureName_b = "fig_n-1_"+target_var+"__"+production_b+"_"+dataType_b+"_"+fileTag_b+"_"+sample_b+"_"+CopyReplaceAll(year_string_b, ",","_");
      if (target_var == "met" && Contains(sample_a, "zll")) {
        pm.Push<Hist1D>(axis_dict[target_var+"_zll"],
          base_filters && n_minus_1_cut_b,
          procs_mc, plt_lin).Weight(weight).Tag("FixName:"+figureName_b);
      } else if (target_var == "met" && Contains(sample_a, "ttbar")) {
        pm.Push<Hist1D>(axis_dict[target_var+"_ttbar"],
          base_filters && n_minus_1_cut_b,
          procs_mc, plt_log).Weight(weight).Tag("FixName:"+figureName_b);
      // For log plots
      } else if (std::find(log_plots.begin(), log_plots.end(), target_var) != log_plots.end()) {
        pm.Push<Hist1D>(axis_dict[target_var],
          base_filters&&n_minus_1_cut_b,
          procs_mc, plt_log).Weight(weight).Tag("FixName:"+figureName_b);
      // Normal case
      } else {
        pm.Push<Hist1D>(axis_dict[target_var],
          base_filters&&n_minus_1_cut_b,
          procs_mc, plt_lin).Weight(weight).Tag("FixName:"+figureName_b);
      }
    }

  }

  //// Example
  //vector<shared_ptr<Process> > procs;
  //NamedFunc sampleCuts_a = "1";
  //for (auto & item : sampleCutDict_a) {
  //  sampleCuts_a = sampleCuts_a && item.value();
  //}
  //NamedFunc sampleCuts_b = "1";
  //for (auto & item : sampleCutDict_b) {
  //  sampleCuts_b = sampleCuts_b && item.value();
  //}
  //pm.Push<Hist1D>(axis_dict["met"],
  //  base_filters,
  //  procs, plt_shapes).Weight(weight).Tag("FixName:fig_n-1_"+sample_a+"_met_"+CopyReplaceAll(year_string_a, ",","_"));

  pm.multithreaded_ = !single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.);

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      //{"sample", required_argument, 0, 0},
      //{"year", required_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      //if(optname == "sample"){
      //  sample = optarg;
      //} else if (optname == "year") {
      //  year_string = optarg;
      if (0) {
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}

