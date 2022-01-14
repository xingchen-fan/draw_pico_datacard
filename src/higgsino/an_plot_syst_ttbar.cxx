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

namespace{
  bool single_thread = false;
  string year_string = "run2";
  bool unblind = true;
  // string_options is split by comma. ex) option1,option2 
  // Use HigUtilities::is_in_string_options(string_options, "option2") to check if in string_options.
  string string_options = "makePies,plot_data_vs_mc,plot_isr,plot_kappa";
  // Other options: plot_planes,plot_btags,plot_isrTrueb,plot_bCorrelation,paper_style
}

void combine_bins(vector<pair<string, NamedFunc> > & combined_bins, vector<pair<string, NamedFunc> > const & bins_a,  vector<pair<string, NamedFunc> > const & bins_b) {
  for (auto const & bin_a : bins_a) {
    for (auto const & bin_b : bins_b) {
      combined_bins.push_back({bin_a.first+"_"+bin_b.first,bin_a.second && bin_b.second});
      //cout<<bin_a.first+"_"+bin_b.first<<" "<<(bin_a.second && bin_b.second).Name()<<endl;
    }
  }
}
void combine_bins(vector<pair<string, NamedFunc> > & combined_bins, vector<pair<string, NamedFunc> > const & bins_a,  vector<pair<string, NamedFunc> > const & bins_b, vector<pair<string, NamedFunc> > const & bins_c) {
  for (auto const & bin_a : bins_a) {
    for (auto const & bin_b : bins_b) {
      for (auto const & bin_c : bins_c) {
        combined_bins.push_back({bin_a.first+"_"+bin_b.first+"_"+bin_c.first, bin_a.second && bin_b.second && bin_c.second});
        //cout<<bin_a.first+"_"+bin_b.first+"_"+bin_c.first<<" "<<(bin_a.second && bin_b.second && bin_c.second).Name()<<endl;
      }
    }
  }
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);
  cout<<"Will plot: "<<string_options<<endl;

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
  PlotOpt lin_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::simulation_supplementary);
  PlotOpt lin_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::supplementary).Bottom(BottomType::ratio);
  PlotOpt log_norm_paper = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::simulation_supplementary);
  PlotOpt log_norm_paper_data = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::supplementary).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_info_paper = lin_norm().Stack(StackType::shapes).Bottom(BottomType::off).Title(TitleType::simulation_supplementary);
  PlotOpt lin_shapes_paper_data = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio).Title(TitleType::supplementary);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_data = {lin_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_log = {log_norm_data};
  if (HigUtilities::is_in_string_options(string_options, "paper_style")) {
     plt_shapes = {lin_shapes_paper};
     plt_shapes_data = {lin_shapes_paper_data};
     if (!unblind) {
       plt_lin = {lin_norm_paper};
       plt_log = {log_norm_paper};
     }
     else {
       plt_lin = {lin_norm_paper_data};
       plt_log = {log_norm_paper_data};
     }
  }

  string higgsino_version = "";

  // Set options
  //string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string mc_base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath_v3/";
  string mc_skim_folder = "mc/merged_higmc_higloose/";
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";

  //string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string data_base_folder = "/net/cms17/cms17r0/pico/NanoAODv7/higgsino_klamath/";
  string data_skim_folder = "data/merged_higdata_higloose/";
  string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";

  //string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string sig_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v3/";
  //string ttbar_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep1T/";
  string ttbar_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed/";

  //years = {2016, 2017, 2018};
  set<int> years;
  HigUtilities::parseYears(year_string, years);
  int lumi_precision = 1;
  if (HigUtilities::is_in_string_options(string_options, "paper_style"))
    lumi_precision = 0;
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string, lumi_precision);

  NamedFunc weight = Higfuncs::final_weight;

  NamedFunc lepton_triggers = Higfuncs::el_trigger || Higfuncs::mu_trigger;
  NamedFunc met_triggers = Higfuncs::met_trigger;;
  NamedFunc triggers_data = lepton_triggers || met_triggers;

  NamedFunc base_filters = Higfuncs::final_ttbar_pass_filters;

  // resolved cuts
  //NamedFunc search_resolved = "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //                       "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //                       "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))"
                         &&Higfuncs::lead_signal_lepton_pt>30;


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
  mctags["other"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  mctags["other_and_single_t"]   = set<string>({"*_WH*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root","*_ST_*.root"});
  // Combine all tags
  mctags["all"] = set<string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  vector<shared_ptr<Process> > ttbar_procs;
  if (HigUtilities::is_in_string_options(string_options, "paper_style")) {
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                      attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other_and_single_t"]),"stitch"));
  }
  else {
    // Set mc processes
    ////ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
    //                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
    //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
    //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["vjets"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["single_t"]),"stitch"));
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other"]),"stitch"));
  }

  if (unblind) {
    ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),triggers_data));
  }

  //// Set signal mc
  //vector<string> signal_mass = {"200", "450", "700", "950"}; 
  //vector<int> sig_colors = {kGreen+1, kRed, kBlue, kOrange}; // need signal_mass.size() >= sig_colors.size()
  //for (unsigned isig(0); isig<signal_mass.size(); isig++){
  //  procs.push_back(Process::MakeShared<Baby_pico>("TChiHH("+signal_mass[isig]+",1)", Process::Type::signal, 
  //    sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+signal_mass[isig]+"_mLSP-0*.root"}), "1"));
  //}

  // Set processes according to btag
  // Getting colors
  // TColor * color;; float red, green, blue;
  // color = gROOT->GetColor(kAzure+1); color->GetRGB(red,green,blue); cout<<red*255<<" "<<green*255<<" "<<blue*255<<endl;
  vector<shared_ptr<Process> > ttbar_procs_btag;
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  ttbar_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));

  // Set processes according to true number of b
  vector<shared_ptr<Process> > ttbar_procs_trueB;
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("0 B-hadron",       Process::Type::background, colors("true_0b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub<1));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("1 B-hadron",       Process::Type::background, colors("true_1b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==1));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("2 B-hadrons",      Process::Type::background, colors("true_2b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub==2));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("3 B-hadrons",      Process::Type::background, colors("true_3b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),  "stitch"&& Functions::ntrub==3));
  ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>
    ("4 or more B-hadrons", Process::Type::background, colors("true_4b"), attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]), "stitch" && Functions::ntrub>=4));
  if (unblind) {
    ttbar_procs_trueB.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    triggers_data
                    ));
  }

  vector<shared_ptr<Process> > ttbar_procs_nisr;
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 0", Process::Type::background,colors("nisr_0"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr==0"));
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 1", Process::Type::background,colors("nisr_1"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr==1"));
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} = 2", Process::Type::background,colors("nisr_2"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr==2"));
  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("N_{ISR} => 3", Process::Type::background,colors("nisr_3"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&nisr>=3"));
  //if (unblind) {
  //  ttbar_procs_nisr.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
  //                  attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
  //                  triggers_data
  //                  ));
  //}

  vector<shared_ptr<Process> > ttbar_data_procs_btag;
  if (unblind) {
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("2b Data", Process::Type::background, colors("0b"),
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt==2&&nbm==2)" && triggers_data
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("3b Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm==3&&nbl==3)" && triggers_data
                    ));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("4b Data", Process::Type::data, kGreen,
                    attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
                    "(nbt>=2&&nbm>=3&&nbl>=4)" && triggers_data
                    ));
  } else {
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (2b)", Process::Type::background,colors("2b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (3b)", Process::Type::background,colors("3b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
    ttbar_data_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. (4b)", Process::Type::background,colors("4b"),
                    attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));
  }

  vector<shared_ptr<Process> > ttbar_mc_procs_btag;
  ttbar_mc_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. MC (2b)", Process::Type::background,colors("2b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt==2&&nbm==2)"));
  ttbar_mc_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. MC (3b)", Process::Type::background,colors("3b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm==3&&nbl==3)"));
  ttbar_mc_procs_btag.push_back(Process::MakeShared<Baby_pico>("All bkg. MC (4b)", Process::Type::background,colors("4b"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["all"]),"stitch&&(nbt>=2&&nbm>=3&&nbl>=4)"));

  PlotMaker pm;

  bool makePies = HigUtilities::is_in_string_options(string_options, "makePies");

  if (makePies) {
    // 2b, (met 0, 75, 150, 200, 300), low drmax
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

    // 3b, (met 0, 75, 150, 200, 300), low drmax
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

    // 4b, (met 0, 75, 150, 200, 300), low drmax
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

    // 2b, (met 0, 75, 150, 200, 300), high drmax
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__2b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

    // 3b, (met 0, 75, 150, 200, 300), high drmax
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__3b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

    // 4b, (met 0, 75, 150, 200, 300), high drmax
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies__4b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs, true, true, true);

    // 2b3b4b
    pm.Push<Table>("FixName:syst__ttbar_pies__2b3b4b"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved, 0, 0, weight)}), ttbar_procs, true, true, true);
  }

  bool plot_data_vs_mc = HigUtilities::is_in_string_options(string_options, "plot_data_vs_mc");

  if (plot_data_vs_mc) {
    // comparison of data with mc
    // [drmax;no drmax] (2b, 3b+4b, all)
    pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "(nbt==2&&nbm==2)"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_drmax_2b").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_drmax_3b4b").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, "hig_cand_drmax[0]", "#DeltaR_{max}", {2.2}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "(nbt>=2)" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_drmax_2b3b4b").LuminosityTag(total_luminosity_string);

    // [<m>;no <m>, low_drmax] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "hig_cand_drmax[0]<1.1&&"
      "(nbt==2&&nbm==2)" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_lowdrmax").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "hig_cand_drmax[0]<1.1&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_lowdrmax").LuminosityTag(total_luminosity_string);
    // [<m>;no <m>, high_drmax] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "hig_cand_drmax[0]>=1.1&&"
      "(nbt==2&&nbm==2)" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_highdrmax").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "hig_cand_drmax[0]>=1.1&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_highdrmax").LuminosityTag(total_luminosity_string);

    // Extra plots [<m>;no <m>] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "(nbt==2&&nbm==2)" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b").LuminosityTag(total_luminosity_string);

    // Supplementary plots [<m>;no <m>;met>75] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "(nbt==2&&nbm==2)&&met>75" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_met75").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&met>75" 
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_met75").LuminosityTag(total_luminosity_string);
    
    // [MET] (all)
    pm.Push<Hist1D>(Axis(16, 0, 800, "met", "p_{T}^{miss} [GeV]", {}),
      base_filters&&ttbar_resolved,
      ttbar_procs, plt_log).Weight(weight).Tag("FixName:syst__ttbar_met").LuminosityTag(total_luminosity_string);


    NamedFunc no_mass_cut = "nlep==1&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<80";
    // [<m>;no <m>, low_drmax, met150-200] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "(nbt==2&&nbm==2)&&" 
      "met>150&&met<=200"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_lowdrmax_met150").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&" 
      "met>150&&met<=200"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_lowdrmax_met150").LuminosityTag(total_luminosity_string);
    // [<m>;no <m>, low_drmax, met200-300] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "(nbt==2&&nbm==2)&&" 
      "met>200&&met<=300",
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_lowdrmax_met200").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&" 
      "met>200&&met<=300"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_lowdrmax_met200").LuminosityTag(total_luminosity_string);
    // [<m>;no <m>, low_drmax, met300-400] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "(nbt==2&&nbm==2)&&" 
      "met>300&&met<=400"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_lowdrmax_met300").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&" 
      "met>300&&met<=400"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_lowdrmax_met300").LuminosityTag(total_luminosity_string);
    // [<m>;no <m>, low_drmax, met400] (2b, 3b+4b)
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "(nbt==2&&nbm==2)&&" 
      "met>400"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_lowdrmax_met400").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&no_mass_cut&&
      "hig_cand_drmax[0]<1.1&&"
      "((nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&" 
      "met>400"
      &&Higfuncs::lead_signal_lepton_pt>30,
      ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b4b_lowdrmax_met400").LuminosityTag(total_luminosity_string);
  }

  bool plot_isr = HigUtilities::is_in_string_options(string_options, "plot_isr");

  if (plot_isr) {
    // [<m> (2b mc, 3b mc, 4b mc); (low, high) drmax]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
      ttbar_mc_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax__mc").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
      ttbar_mc_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax__mc").LuminosityTag(total_luminosity_string);

    // [<m> (2b data, 3b data, 4b data); (low, high) drmax]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
      ttbar_data_procs_btag, plt_shapes_data).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax__data").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
      ttbar_data_procs_btag, plt_shapes_data).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax__data").LuminosityTag(total_luminosity_string);

    // For mc shapes
    // [nisr shapes (2b,3b,4b);(low, high) drmax]
    pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr", "N_{ISR jets}", {}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_nisr__lowdrmax").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr", "N_{ISR jets}", {}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_nisr__highdrmax").LuminosityTag(total_luminosity_string);
    // [<m> shapes (nisr=0,1,2);(low, high) drmax]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
      ttbar_procs_nisr, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_nisr__lowdrmax").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
      ttbar_procs_nisr, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_nisr__highdrmax").LuminosityTag(total_luminosity_string);
    // [<m> shapes (2b, 3b, 4b);(low, high) drmax]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax").LuminosityTag(total_luminosity_string);
    // [<m> shapes (2b, 3b, 4b);very high drmax]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&
      "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
      "hig_cand_drmax[0]>=2.2",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__veryhighdrmax").LuminosityTag(total_luminosity_string);
    // [<m> shapes (2b, 3b, 4b);low drmax, met (0,75,150,200,300)]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>0&&met<=75",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_0").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>75&&met<=150",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_75").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>150&&met<=200",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_150").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>200&&met<=300",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_200").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&met>300",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_met_300").LuminosityTag(total_luminosity_string);
    // [<m> shapes (2b, 3b, 4b);high drmax, met (0,75,150,200,300)]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>0&&met<=75",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_0").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>75&&met<=150",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_75").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>150&&met<=200",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_150").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>200&&met<=300",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_200").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&met>300",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_met_300").LuminosityTag(total_luminosity_string);
    // [<m> shapes (2b, 3b, 4b);low drmax, nisr (0,1,2)]
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&nisr==0",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_nisr_0").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&nisr==1",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_nisr_1").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      base_filters&&ttbar_resolved&&"hig_cand_drmax[0]<1.1&&nisr==2",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__lowdrmax_nisr_2").LuminosityTag(total_luminosity_string);
    // [<m> shapes (2b, 3b, 4b);high drmax, nisr (0,1,2)]
    // [TODO] Strange that drmax[1.1,2.2] shows a correlation. Doesn't agree well with nisr explanation.
    // Could study with xy, xy. Is b[pt high=x, high low=y] close to t1b, t1w, t2b, t2b, none==isr? => bb = 16 combinations, bbbb = 256 combinations
    // Could study with xy, xy. Is b[pt high=x, high low=y] close to t1(b,w), t2(b,w), none==isr? => bb = 9 combinations, bbbb = 81 combinations
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      //base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&nisr==0",
      base_filters&&
      "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
      "hig_cand_drmax[0]>=1.1&&nisr==0",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_nisr_0").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      //base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&nisr==1",
      base_filters&&
      "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
      "hig_cand_drmax[0]>=1.1&&nisr==1",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_nisr_1").LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
      //base_filters&&ttbar_resolved&&"hig_cand_drmax[0]>=1.1&&nisr==2",
      base_filters&&
      "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
      "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
      "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))&&"
      "hig_cand_drmax[0]>=1.1&&nisr==2",
      ttbar_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__highdrmax_nisr_2").LuminosityTag(total_luminosity_string);
  }

  bool plot_kappa = HigUtilities::is_in_string_options(string_options, "plot_kappa");
  if (plot_kappa) {
    if (unblind) {
      system(("./run/higgsino/plot_kappas.exe --sample ttbar --scen data --unblind --year "+year_string).c_str());
    } else {
      system(("./run/higgsino/plot_kappas.exe --sample ttbar --scen mc_as_data --year "+year_string).c_str());
    }
  }

  bool plot_planes = HigUtilities::is_in_string_options(string_options, "plot_planes");
  if (plot_planes && unblind) {
    torch::OrderedDict<string, NamedFunc> regions;
    regions.insert("lowdrmax_met0", "hig_cand_drmax[0]<=1.1&&met>0&&met<=75");
    regions.insert("lowdrmax_met75", "hig_cand_drmax[0]<=1.1&&met>75&&met<=150");
    regions.insert("lowdrmax_met150","hig_cand_drmax[0]<=1.1&&met>150&&met<=200");
    regions.insert("lowdrmax_met200","hig_cand_drmax[0]<=1.1&&met>200&&met<=300");
    regions.insert("lowdrmax_met300","hig_cand_drmax[0]<=1.1&&met>300");
    regions.insert("lowdrmax_met0above","hig_cand_drmax[0]<=1.1&&met>0");
    regions.insert("lowdrmax_met150above","hig_cand_drmax[0]<=1.1&&met>150");
    regions.insert("highdrmax_met0", "hig_cand_drmax[0]>1.1&&met>0&&met<=75");
    regions.insert("highdrmax_met75", "hig_cand_drmax[0]>1.1&&met>75&&met<=150");
    regions.insert("highdrmax_met150","hig_cand_drmax[0]>1.1&&met>150&&met<=200");
    regions.insert("highdrmax_met200","hig_cand_drmax[0]>1.1&&met>200&&met<=300");
    regions.insert("highdrmax_met300","hig_cand_drmax[0]>1.1&&met>300");
    regions.insert("highdrmax_met0above","hig_cand_drmax[0]>1.1&&met>0");
    regions.insert("highdrmax_met150above","hig_cand_drmax[0]>1.1&&met>150");
    for (auto const & region : regions) {
      // [<m> (2b data, 3b data, 4b data); regions]
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        base_filters&&ttbar_resolved&&region.value(),
        ttbar_data_procs_btag, plt_shapes).Weight(weight).Tag("FixName:syst__ttbar_amjj_btag__region_"+region.key()+"__data").LuminosityTag(total_luminosity_string);
      // [<m>;no <m>, low_drmax, met400] (2b, 3b, 4b)
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)"&&region.value(),
        ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_2b_region_"+region.key()).LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)"&&region.value(),
        ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_3b_region_"+region.key()).LuminosityTag(total_luminosity_string);
      pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)"&&region.value(),
        ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_4b_region_"+region.key()).LuminosityTag(total_luminosity_string);
    }
  }

  bool plot_btags = HigUtilities::is_in_string_options(string_options, "plot_btags");
  if (plot_btags && unblind) {
    // Plot in bins
    vector<pair<string, NamedFunc> > mass_bins = {{"SDB","hig_cand_am[0]<=100||hig_cand_am[0]>140"},{"HIG","hig_cand_am[0]>100&&hig_cand_am[0]<=140"}};
    vector<pair<string, NamedFunc> > met_bins = {{"MET75", "met<=75"},{"75MET150","met>75&&met<=150"},{"150MET200","met>150&&met<=200"},{"200MET", "met>200"}};
    vector<pair<string, NamedFunc> > drmax_bins = {{"low-drmax","hig_cand_drmax[0]<=1.1"}, {"high-drmax","hig_cand_drmax[0]>1.1"}};
    // combine bins
    vector<pair<string, NamedFunc> > met_drmax_mass_bins;
    combine_bins(met_drmax_mass_bins, drmax_bins, met_bins, mass_bins);
    // Plot nb in each bin
    for (auto const & bin : met_drmax_mass_bins) {
      pm.Push<Hist1D>(Axis(3, 1.5, 4.5, Higfuncs::hig_bcat, "N_{b}", {2.5}),
        base_filters&&ttbar_resolved&&bin.second, 
        ttbar_procs_trueB, plt_log).Weight(weight).Tag("FixName:syst__ttbar_btags_"+bin.first+"_"+CopyReplaceAll(year_string, ",","_")).LuminosityTag(total_luminosity_string);
    }

    torch::OrderedDict<string, NamedFunc> drmaxs;
    drmaxs.insert("lowdrmax", "hig_cand_drmax[0]<=1.1");
    drmaxs.insert("highdrmax", "hig_cand_drmax[0]>1.1");
    // (all, 2b, 3b, 4b)
    torch::OrderedDict<string, NamedFunc> btags;
    btags.insert("2b", "nbt==2&&nbm==2");
    btags.insert("3b", "nbt>=2&&nbm==3&&nbl==3");
    btags.insert("4b", "nbt>=2&&nbm>=3&&nbl>=4");
    for (auto const & drmax : drmaxs) {
      for (auto const & btag : btags) {
        // <m_bb> across true 0b, 1b, 2b, 3b, 4b for met>150
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          base_filters&&ttbar_resolved&&drmax.value()&&btag.value()&&"met>150",
          ttbar_procs_trueB, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_trueb01234_"+drmax.key()+"_"+btag.key()+"_met150").LuminosityTag(total_luminosity_string);
      }
    }

  }

  bool plot_isrTrueb = HigUtilities::is_in_string_options(string_options, "plot_isrTrueb");
  if (plot_isrTrueb) {
    torch::OrderedDict<string, NamedFunc> drmaxs;
    drmaxs.insert("lowdrmax", "hig_cand_drmax[0]<=1.1");
    drmaxs.insert("highdrmax", "hig_cand_drmax[0]>1.1");

    // (all, 2b, 3b, 4b)
    torch::OrderedDict<string, NamedFunc> btags;
    btags.insert("2b3b4b", "(nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)");
    btags.insert("2b", "nbt==2&&nbm==2");
    btags.insert("3b", "nbt>=2&&nbm==3&&nbl==3");
    btags.insert("4b", "nbt>=2&&nbm>=3&&nbl>=4");

    for (auto const & drmax : drmaxs) {
      for (auto const & btag : btags) {
        // <m_bb> across ISR=0,ISR=1,ISR=2
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          base_filters&&ttbar_resolved&&drmax.value()&&btag.value(),
          ttbar_procs_nisr, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_isr012_"+drmax.key()+"_"+btag.key()).LuminosityTag(total_luminosity_string);
        // <m_bb> across ISR=0,ISR=1,ISR=2 for met>150
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          base_filters&&ttbar_resolved&&drmax.value()&&btag.value()&&"met>150",
          ttbar_procs_nisr, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_isr012_"+drmax.key()+"_"+btag.key()+"_met150").LuminosityTag(total_luminosity_string);
        // <m_bb> across true 0b, 1b, 2b, 3b, 4b
        pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
          base_filters&&ttbar_resolved&&drmax.value()&&btag.value(),
          ttbar_procs_trueB, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_trueb01234_"+drmax.key()+"_"+btag.key()).LuminosityTag(total_luminosity_string);
        //// <m_bb> across true 0b, 1b, 2b, 3b, 4b for met>150
        //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
        //  base_filters&&ttbar_resolved&&drmax.value()&&btag.value()&&"met>150",
        //  ttbar_procs_trueB, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_amjj_trueb01234_"+drmax.key()+"_"+btag.key()+"_met150").LuminosityTag(total_luminosity_string);
      }
    }
  }

  bool plot_bCorrelation = HigUtilities::is_in_string_options(string_options, "plot_bCorrelation");
  if (plot_bCorrelation) {
    // (all, 2b, 3b, 4b)
    torch::OrderedDict<string, NamedFunc> btags;
    btags.insert("2b3b4b", "(nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4)");
    btags.insert("2b", "nbt==2&&nbm==2");
    btags.insert("3b", "nbt>=2&&nbm==3&&nbl==3");
    btags.insert("4b", "nbt>=2&&nbm>=3&&nbl>=4");
    NamedFunc base_cut = "met>0";
    for (auto const & btag : btags) {
      // pt of jets
      pm.Push<Hist1D>(Axis(16, 0, 480., "jet_pt", "p_{T} jet [GeV]", {}),
        base_filters&&ttbar_resolved&&base_cut&&btag.value(),
        ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_jet_pt_"+btag.key()).LuminosityTag(total_luminosity_string);
      // dr of h1
      pm.Push<Hist1D>(Axis(20,0,4,Higfuncs::h1_dr, "Lead higgs #DeltaR", {1.1, 2.2}),
        base_filters&&ttbar_resolved&&base_cut&&btag.value(),
        ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_h1_dr_"+btag.key()).LuminosityTag(total_luminosity_string);
      // dr of h2
      pm.Push<Hist1D>(Axis(20,0,4,Higfuncs::h2_dr, "Sublead higgs #DeltaR", {1.1, 2.2}),
        base_filters&&ttbar_resolved&&base_cut&&btag.value(),
        ttbar_procs, plt_lin).Weight(weight).Tag("FixName:syst__ttbar_h2_dr_"+btag.key()).LuminosityTag(total_luminosity_string);
    }
  }

  if (makePies) {
    // According to true b
    // 2b, (met 0, 75, 150, 200, 300), low drmax 
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

    // 3b, (met 0, 75, 150, 200, 300), low drmax
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

    // 4b, (met 0, 75, 150, 200, 300), low drmax
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met0_lowdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met75_lowdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met150_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met200_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met300_lowdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

    // 2b, (met 0, 75, 150, 200, 300), high drmax
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__2b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt==2&&nbm==2)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

    // 3b, (met 0, 75, 150, 200, 300), high drmax
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__3b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);

    // 4b, (met 0, 75, 150, 200, 300), high drmax
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met0_highdrmax"  , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)         &&met<=75 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met75_highdrmax" , vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>75 &&met<=150&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met150_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150&&met<=200&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met200_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200&&met<=300&&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
    pm.Push<Table>("FixName:syst__ttbar_pies_trueb__4b_met300_highdrmax", vector<TableRow> ({TableRow("", base_filters&&ttbar_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300          &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), ttbar_procs_trueB, true, true, true);
  }

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
      {"unblind", no_argument, 0, 'u'},
      {"year", required_argument, 0, 'y'},
      {"string_options", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "so:uy:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 'u':
      unblind = true;
      break;
    case 'y':
      year_string = optarg;
      break;
    case 'o':
      string_options = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
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
