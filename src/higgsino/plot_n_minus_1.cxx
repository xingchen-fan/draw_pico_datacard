#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "higgsino/ordered_dict.hpp"

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
#include "core/hist2d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  bool unblind = false;
  // sample can be search, ttbar, zll, qcd
  string sample_name = "search";
  // year can be 2016, 2017, 2018 or 2016,2017,2018 or run2
  string year_string = "2016";
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

int main(int argc, char *argv[]){
  // ./run/higgsino/plot_n_minus_1.exe --sample search --year 2016
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2);
  PlotOpt log_norm_data = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.2).Bottom(BottomType::ratio);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> plt_2D = {style().Stack(StackType::data_norm).Title(TitleType::data)};

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_log = {log_norm_data};

  // Set options
  string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  //string mc_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado";
  // preselect:
  //   ((nbt>=2 && njet>=4 && njet<=5)||(Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1))
  //   nvlep==0 && ntk==0 && !low_dphi_met && met>150 && 
  // higloose: 
  //   (nbt>=2 || nbdft>=2 || Sum$(fjet_pt>300 && fjet_msoftdrop>50)>0)&&
  //   met>150 && nvlep==0
  string search_mc_skim_folder = "mc/merged_higmc_higloose/";
  // higlep1T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || ((nbt>=2 || nbdft>=2) && njet>=4 && njet<=5)) &&
  //   nlep==1 && 
  //   (Max$(el_pt*el_sig)>40 || Max$(mu_pt*mu_sig)>40) // pass_1l_trig40 (sig is signal lepton)
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  // higlep2T:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nlep==2 && 
  //   @ll_m.size()>=1 && Sum$(ll_m>80 && ll_m<100)>=1
  //   (Max$(el_pt*el_sig)>30 || Max$(mu_pt*mu_sig)>30) // pass_2l_trig30
  string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  // higqcd with met150:
  //   (Sum$(fjet_pt>300 && fjet_msoftdrop>50)>1 || (njet>=4 && njet<=5))
  //   nvlep==0 && ntk==0 && low_dphi_met &&
  //   met>150  // Since applied to met150 skim
  string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt";
  string search_data_skim_folder = "data/merged_higmc_higloose/";
  string ttbar_data_skim_folder = "data/merged_higmc_higlep1T/";
  string zll_data_skim_folder = "data/merged_higmc_higlep2T/";
  string qcd_data_skim_folder = "data/merged_higmc_higqcd/";

  //string sig_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_preselect/";

  set<int> years;
  HigUtilities::parseYears(year_string, years);
  //years = {2016, 2017, 2018};
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig*w_years;
  //if (years.size()==1 && *years.begin()==2016) weight *= "137.";
  //else weight *= w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig*w_years;
  //NamedFunc weight = "weight"*Higfuncs::eff_higtrig_run2*w_years;

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
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

  string mc_skim_folder;
  if (sample_name == "ttbar") mc_skim_folder = ttbar_mc_skim_folder;
  else if (sample_name == "zll") mc_skim_folder = zll_mc_skim_folder;
  else if (sample_name == "qcd") mc_skim_folder = qcd_mc_skim_folder;
  else mc_skim_folder = search_mc_skim_folder;
  string data_skim_folder;
  if (sample_name == "ttbar") data_skim_folder = ttbar_data_skim_folder;
  else if (sample_name == "zll") data_skim_folder = zll_data_skim_folder;
  else if (sample_name == "qcd") data_skim_folder = qcd_data_skim_folder;
  else data_skim_folder = search_data_skim_folder;

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
  torch::OrderedDict<string, NamedFunc> ttbar_resolved_cuts;
  ttbar_resolved_cuts.insert("nlep", "nlep==1");
  //ttbar_resolved_cuts.insert("lep_pt", "lep_pt[0]>30");
  ttbar_resolved_cuts.insert("lep_pt", Higfuncs::lead_signal_lepton_pt>30);
  ttbar_resolved_cuts.insert("mt", "mt<=100");
  ttbar_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  ttbar_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  ttbar_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  ttbar_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  ttbar_resolved_cuts.insert("btags", "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))");
  torch::OrderedDict<string, NamedFunc> zll_resolved_cuts;
  zll_resolved_cuts.insert("nlep", "nlep==2");
  zll_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  zll_resolved_cuts.insert("met", "met<50");
  zll_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  zll_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  zll_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  zll_resolved_cuts.insert("dbtags", "(nbm==0||nbm==1||nbm==2||nbm>=3)");
  torch::OrderedDict<string, NamedFunc> qcd_resolved_cuts;
  qcd_resolved_cuts.insert("low_dphi_met", "low_dphi_met");
  qcd_resolved_cuts.insert("nvlep", "nvlep==0");
  qcd_resolved_cuts.insert("met", "met>150");
  qcd_resolved_cuts.insert("njet", "njet>=4&&njet<=5");
  qcd_resolved_cuts.insert("hig_cand_drmax", "hig_cand_drmax[0]<2.2");
  qcd_resolved_cuts.insert("hig_cand_am", "hig_cand_am[0]<200");
  qcd_resolved_cuts.insert("hig_cand_dm", "hig_cand_dm[0]<40");
  qcd_resolved_cuts.insert("dbtags", "(nbm==0||nbm==1||nbm==2||nbm>=3)");

  vector<shared_ptr<Process> > procs;
  // Set mc processes
  //procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  vector<string> sigm = {"175", "500", "950"};
  //vector<string> sigm = {"400"};
  //vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  vector<int> sig_colors = {kSpring, kPink, kAzure, kOrange}; // need sigm.size() >= sig_colors.size()
  for (unsigned isig(0); isig<sigm.size(); isig++){
    procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
      sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+sigm[isig]+"_mLSP-0*.root"}), 
      "stitch"));
  }

  if (unblind) {
    procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, data_skim_folder, {"*.root"}),"stitch"));
  }

  NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2" && "nlep==1"; //since pass_fsjets is not quite usable...
  //NamedFunc base_filters = Functions::hem_veto && "pass && met/mht<2 && met/met_calo<2";//HigUtilities::pass_2016; //since pass_fsjets is not quite usable...

  PlotMaker pm;

  torch::OrderedDict<string, Axis> axis_dict;
  axis_dict.insert("ntk",Axis(10, -0.5, 9.5, "ntk", "Number of iso tk", {0.5}));
  axis_dict.insert("low_dphi_met",Axis(2, -0.5, 1.5, "low_dphi_met", "Low #Delta #phi", {0.5}));
  axis_dict.insert("nvlep",Axis(6, -0.5, 5.5, "nvlep", "N_{veto leps}", {0.5}));
  axis_dict.insert("nlep",Axis(5, 0.5, 5.5, "nlep", "N_{leps}", {}));
  axis_dict.insert("met",Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}));
  axis_dict.insert("ht",Axis(20, 0, 1000., "ht", "HT [GeV]", {}));
  axis_dict.insert("h1b1_pt",Axis(16, 0, 480., h1b1_pt, "p_{T} lead higgs high-tag b jet [GeV]", {}));
  axis_dict.insert("h1b2_pt",Axis(16, 0, 480., h1b2_pt, "p_{T} lead higgs low-tag b jet [GeV]", {}));
  axis_dict.insert("h2b1_pt",Axis(16, 0, 480., h2b1_pt, "p_{T} sublead higgs high-tag b jet [GeV]", {}));
  axis_dict.insert("h2b2_pt",Axis(16, 0, 480., h2b2_pt, "p_{T} sublead higgs low-tag b jet [GeV]", {}));
  axis_dict.insert("h1_dr",Axis(20,0,4,h1_dr, "Lead higgs #DeltaR", {1.1, 2.2}));
  axis_dict.insert("h2_dr",Axis(20,0,4,h2_dr, "Sublead higgs #DeltaR", {1.1, 2.2}));
  axis_dict.insert("h1_mass",Axis(10, 0, 200, h1_mass, "Lead m_{bb} [GeV]", {100, 140}));
  axis_dict.insert("h2_mass",Axis(10, 0, 200, h2_mass, "Sublead m_{bb} [GeV]", {100, 140}));
  axis_dict.insert("jet_pt[0]",Axis(16, 0, 480., "jet_pt[0]", "p_{T} jet [0] [GeV]", {}));
  axis_dict.insert("jet_pt[1]",Axis(16, 0, 480., "jet_pt[1]", "p_{T} jet [1] [GeV]", {}));
  axis_dict.insert("jet_pt[2]",Axis(16, 0, 480., "jet_pt[2]", "p_{T} jet [2] [GeV]", {}));
  axis_dict.insert("jet_pt[3]",Axis(16, 0, 480., "jet_pt[3]", "p_{T} jet [3] [GeV]", {}));
  axis_dict.insert("jet_pt[4]",Axis(16, 0, 480., "jet_pt[4]", "p_{T} jet [4] [GeV]", {}));
  axis_dict.insert("met_zll",Axis(20, 0, 150., "met", "p_{T}^{miss} [GeV]", {30}));
  axis_dict.insert("njet",Axis(12, -0.5, 11.5, "njet", "N_{jets}", {3.5, 5.5}));
  axis_dict.insert("hig_cand_drmax",Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}));
  axis_dict.insert("hig_cand_am",Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}));
  axis_dict.insert("hig_cand_dm",Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}));
  axis_dict.insert("btags",Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}));
  //axis_dict.insert("lep_pt",Axis(10, 0, 300., "lep_pt[0]", "p_{l} [GeV]", {30}));
  axis_dict.insert("lep_pt",Axis(10, 0, 300., Higfuncs::lead_signal_lepton_pt, "p_{l} [GeV]", {30}));
  axis_dict.insert("mt",Axis(10, 0, 200., "mt", "m_{T} [GeV]", {100}));
  axis_dict.insert("dbtags",Axis(5, -0.5, 4.5, "nbm", "N_{b medium}", {}));
  axis_dict.insert("ll_pt",Axis(16, 0, 400., "ll_pt[0]", "p_{ll} [GeV]", {75., 150., 200., 300}));
  vector<string> log_plots = {"met", "btags"};


  // Draw n-1
  // Selection according to sample
  torch::OrderedDict<string, NamedFunc> map_resolved_cuts;
  if (sample_name == "ttbar") map_resolved_cuts = ttbar_resolved_cuts;
  else if (sample_name == "zll") map_resolved_cuts = zll_resolved_cuts;
  else if (sample_name == "qcd") map_resolved_cuts = qcd_resolved_cuts;
  else map_resolved_cuts = search_resolved_cuts;

  // Variables to draw
  set<string> target_variables;
  // Add all cut variables
  for (auto const & target_var : map_resolved_cuts) {
    target_variables.insert(target_var.key());
  }
  // Add additional variables
  target_variables.insert("ht");
  target_variables.insert("jet_pt[0]");
  target_variables.insert("jet_pt[1]");
  target_variables.insert("jet_pt[2]");
  target_variables.insert("jet_pt[3]");
  target_variables.insert("jet_pt[4]");
  target_variables.insert("h1b1_pt");
  target_variables.insert("h1b2_pt");
  target_variables.insert("h2b1_pt");
  target_variables.insert("h2b2_pt");
  target_variables.insert("h1_dr");
  target_variables.insert("h2_dr");
  target_variables.insert("h1_mass");
  target_variables.insert("h2_mass");
  if (sample_name=="zll") target_variables.insert("ll_pt");
  if (sample_name=="ttbar") target_variables.insert("met");


  // Loop over variables
  for (auto const & target_var : target_variables) {
    // Make n-1 cut for target_var
    NamedFunc n_minus_1_cut = "1";
    for (auto const & cut_item : map_resolved_cuts) {
      // Remove plotting variable
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
      n_minus_1_cut = n_minus_1_cut && cut_item.value();
    }
    
    // Draw target_var
    // For special case
    if (target_var == "met" && sample_name == "zll") {
      pm.Push<Hist1D>(axis_dict[target_var+"_zll"],
        base_filters&&n_minus_1_cut,
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    } else if (target_var == "jet_pt[4]") {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut&&"njet>=5",
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }
    // For log plots
    else if (std::find(log_plots.begin(), log_plots.end(), target_var) != log_plots.end()) {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut,
        procs, plt_log).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    // Normal case
    } else {
      pm.Push<Hist1D>(axis_dict[target_var],
        base_filters&&n_minus_1_cut,
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_"+target_var+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }
  }

  NamedFunc resolved_cuts = "1";
  for (auto & item : map_resolved_cuts) {
    resolved_cuts = resolved_cuts && item.value();
  }
  // Draw 2D ht vs met
  pm.Push<Hist2D>(axis_dict["ht"], axis_dict["met"],
    base_filters&&resolved_cuts, procs, plt_2D).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ht_vs_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);

  //// Special case for zll. Draw but do not cut on ll_pt
  //if (sample_name=="zll") {
  //  pm.Push<Hist1D>(Axis(16, 0, 400., "ll_pt[0]", "p_{ll} [GeV]", {75., 150., 200., 300}),
  //    base_filters&&resolved_cuts,
  //    procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_ll_pt_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  //}
  //// Special case for ttbar. Draw but do not cut on met
  //if (sample_name=="ttbar") {
  //  pm.Push<Hist1D>(axis_dict["met"],
  //    base_filters&&resolved_cuts,
  //    procs, plt_log).Weight(weight).Tag("FixName:fig_n-1_"+sample_name+"_met_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
  //}

  //// Example
  //NamedFunc resolved_cuts = "1";
  //torch::OrderedDict<string, NamedFunc> map_resolved_cuts;
  //if (sample_name == "ttbar") map_resolved_cuts = ttbar_resolved_cuts;
  //else if (sample_name == "zll") map_resolved_cuts = zll_resolved_cuts;
  //else if (sample_name == "qcd") map_resolved_cuts = qcd_resolved_cuts;
  //else map_resolved_cuts = search_resolved_cuts;
  //for (auto & item : map_resolved_cuts) {
  //  resolved_cuts = resolved_cuts && item.value();
  //}
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters&&resolved_cuts&&"hig_cand_drmax[0]<1.1",
  //  procs, plt_lin).Weight(weight).Tag("FixName:n-1_hig_cand_am_2016");

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
      {"unblind", no_argument, 0, 0},
      {"sample", required_argument, 0, 0},
      {"year", required_argument, 0, 0},
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
      if(optname == "sample"){
        sample_name = optarg;
      } else if (optname == "year") {
        year_string = optarg;
      } else if (optname == "unblind") {
        unblind = true;
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

