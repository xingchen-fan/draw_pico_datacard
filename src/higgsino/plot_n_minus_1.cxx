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
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

const NamedFunc w_years("w_years", [](const Baby &b) -> NamedFunc::ScalarType{
  if (b.SampleType()<0) return 1.;

  double weight = 1;
  if (b.type()==106000) {
    return 35.9;
  }
  if (b.SampleType()==2016){
    return weight*35.9;
  } else if (b.SampleType()==2017){
    return weight*41.5;
  } else {
    return weight*59.6;
  }
});

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
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
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info);
  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};

  // Set options
  //string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt/";
  string mc_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado";
  string search_mc_skim_folder = "mc/merged_higmc_higloose/";
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  string data_base_folder = "/net/cms25/cms25r0/pico/NanoAODv5/higgsino_humboldt";
  string search_data_skim_folder = "data/merged_higmc_higloose/";
  string ttbar_data_skim_folder = "data/merged_higmc_higlep1T/";
  string zll_data_skim_folder = "data/merged_higmc_higlep2T/";
  string qcd_data_skim_folder = "data/merged_higmc_higqcd/";

  string sig_base_folder = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_eldorado/";
  string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";

  set<int> years;
  //years = {2016, 2017, 2018};
  years = {2016};

  NamedFunc weight = "w_lumi*w_isr"*Higfuncs::eff_higtrig;
  if (years.size()==1 && *years.begin()==2016) weight *= "137.";
  else weight *= w_years;

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

  // sample can be search, ttbar, zll, qcd
  string sample = "qcd";
  string year_string = "2016";


  string mc_skim_folder;
  if (sample == "ttbar") mc_skim_folder = ttbar_mc_skim_folder;
  else if (sample == "zll") mc_skim_folder = zll_mc_skim_folder;
  else if (sample == "qcd") mc_skim_folder = qcd_mc_skim_folder;
  else mc_skim_folder = search_mc_skim_folder;


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
  ttbar_resolved_cuts.insert("lep_pt", "lep_pt[0]>30");
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
  //procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),"stitch"));


  // Set processes according to btag
  // Getting colors
  // TColor * color; float red, green, blue;
  // color = gROOT->GetColor(kAzure+1); color->GetRGB(red,green,blue); cout<<red*255<<" "<<green*255<<" "<<blue*255<<endl;
  vector<shared_ptr<Process> > procs_btag;
  procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
  procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbm==2)"));

  // Set processes according to true number of b
  vector<shared_ptr<Process> > procs_trueB;
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("3 B-hadrons",      Process::Type::background, colors("true_3b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),  "stitch"&& Functions::ntrub==3));
  procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("4 B-hadrons", Process::Type::background, colors("true_4b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==4));

  // Set processes according to true number of b simple version
  vector<shared_ptr<Process> > procs_trueB012;
  procs_trueB012.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  procs_trueB012.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  procs_trueB012.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));

  vector<shared_ptr<Process> > qcd_data_procs_btag;
  qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (0b)", Process::Type::background,colors("0b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbm==0)"));
  qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (1b)", Process::Type::background,colors("1b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbm==1)"));
  qcd_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&(nbm==2)"));

  vector<shared_ptr<Process> > procs_nisr;
  procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 0", Process::Type::background,colors("nisr_0"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&nisr==0"));
  procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 1", Process::Type::background,colors("nisr_1"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&nisr==1"));
  procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 2", Process::Type::background,colors("nisr_2"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["all"]),"stitch&&nisr==2"));

  NamedFunc base_filters = HigUtilities::pass_2016 && "met/mht<2 && met/met_calo<2"; //since pass_fsjets is not quite usable...

  PlotMaker pm;

  torch::OrderedDict<string, Axis> axis_dict;
  axis_dict.insert("ntk",Axis(10, -0.5, 9.5, "ntk", "Number of iso tk", {0.5}));
  axis_dict.insert("low_dphi_met",Axis(2, -0.5, 1.5, "low_dphi_met", "Low #Delta #phi", {0.5}));
  axis_dict.insert("nvlep",Axis(6, -0.5, 5.5, "nvlep", "N_{veto leps}", {0.5}));
  axis_dict.insert("nlep",Axis(5, 0.5, 5.5, "nvlep", "N_{leps}", {1.5}));
  axis_dict.insert("met",Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}));
  axis_dict.insert("met_zll",Axis(20, 0, 150., "met", "p_{T}^{miss} [GeV]", {30}));
  axis_dict.insert("njet",Axis(12, -0.5, 11.5, "njet", "N_{jets}", {3.5, 5.5}));
  axis_dict.insert("hig_cand_drmax",Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}));
  axis_dict.insert("hig_cand_am",Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}));
  axis_dict.insert("hig_cand_dm",Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}));
  axis_dict.insert("btags",Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}));
  axis_dict.insert("lep_pt",Axis(10, 0, 300., "lep_pt[0]", "p_{l} [GeV]", {30}));
  axis_dict.insert("mt",Axis(10, 0, 200., "mt", "m_{T} [GeV]", {100}));
  axis_dict.insert("dbtags",Axis(5, -0.5, 4.5, "nbm", "N_{b medium}", {}));
  //axis_dict.insert("ll_pt",Axis(16, 0, 400., "ll_pt[0]", "p_{ll} [GeV]", {75., 150., 200., 300}));


  // Draw n-1
  // Selection according to sample
  torch::OrderedDict<string, NamedFunc> map_resolved_cuts;
  if (sample == "ttbar") map_resolved_cuts = ttbar_resolved_cuts;
  else if (sample == "zll") map_resolved_cuts = zll_resolved_cuts;
  else if (sample == "qcd") map_resolved_cuts = qcd_resolved_cuts;
  else map_resolved_cuts = search_resolved_cuts;
  // Loop over variables
  for (auto const & target_var : map_resolved_cuts) {
    // Make n-1 cut for target_var
    NamedFunc n_minus_1_cut = "1";
    for (auto const & cut_item : map_resolved_cuts) {
      if (cut_item.key() == target_var.key()) continue;
      // Below are special conditions
      if (target_var.key() == "njet") {
        if (cut_item.key() == "hig_cand_drmax") continue;
        if (cut_item.key() == "hig_cand_am") continue;
        if (cut_item.key() == "hig_cand_dm") continue;
      }
      if (target_var.key() == "nlep") {
        if (cut_item.key() == "mt") continue;
      }
      n_minus_1_cut = n_minus_1_cut && cut_item.value();
    }
    
    // Draw target_var
    // For special case
    if (target_var.key() == "met" && sample == "zll") {
      pm.Push<Hist1D>(axis_dict[target_var.key()+"_zll"],
        base_filters&&n_minus_1_cut,
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+target_var.key()+"_"+year_string);
    }
    // For log plots
    else if (target_var.key() == "met" || target_var.key() == "btags") {
      pm.Push<Hist1D>(axis_dict[target_var.key()],
        base_filters&&n_minus_1_cut,
        procs, plt_log).Weight(weight).Tag("FixName:fig_n-1_"+target_var.key()+"_"+year_string);
    // Normal case
    } else {
      pm.Push<Hist1D>(axis_dict[target_var.key()],
        base_filters&&n_minus_1_cut,
        procs, plt_lin).Weight(weight).Tag("FixName:fig_n-1_"+target_var.key()+"_"+year_string);
    }
  }

  NamedFunc resolved_cuts = "1";
  for (auto & item : map_resolved_cuts) {
    resolved_cuts = resolved_cuts && item.value();
  }
  // Special case for zll. Draw but do not cut on ll_pt
  if (sample=="zll") {
    pm.Push<Hist1D>(Axis(16, 0, 400., "ll_pt[0]", "p_{ll} [GeV]", {75., 150., 200., 300}),
      base_filters&&resolved_cuts,
      procs, plt_lin).Weight(weight).Tag("FixName:n-1_hig_cand_am_2016");
  }


  //// Example
  //NamedFunc resolved_cuts = "1";
  //torch::OrderedDict<string, NamedFunc> map_resolved_cuts;
  //if (sample == "ttbar") map_resolved_cuts = ttbar_resolved_cuts;
  //else if (sample == "zll") map_resolved_cuts = zll_resolved_cuts;
  //else if (sample == "qcd") map_resolved_cuts = qcd_resolved_cuts;
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
      if(false){
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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
