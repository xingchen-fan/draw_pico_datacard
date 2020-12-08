//this script generates the plots seen in the kinematic variables and resolved selection section of the AN

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
#include "higgsino/apply_trigeffs2016.hpp"
#include "higgsino/apply_trigeffs2017.hpp"
#include "higgsino/apply_trigeffs2018.hpp"
#include "core/cross_sections.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

using namespace std;
using namespace Higfuncs;

namespace{
  bool single_thread = false;
  string year_string = "2016";
  bool do_kinematic_plots = true;
  bool do_nminusone_plots = true;
  bool do_cutflow = true;
  bool do_piecharts = true;
  bool do_srambb_plots = true;
  bool unblind = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  GetOptions(argc, argv);

  Palette colors("txt/colors.txt", "default");

  //------------------------------------------------------------------------------------
  //                                 plot opts
  //------------------------------------------------------------------------------------

  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(PlotOptTypes::TitleType::info)   
    .Bottom(PlotOptTypes::BottomType::off)
    .YAxis(PlotOptTypes::YAxisType::linear)
    .Stack(PlotOptTypes::StackType::signal_overlay).LegendColumns(3);
  PlotOpt log_norm_info = lin_norm_info().YAxis(PlotOptTypes::YAxisType::log);
  PlotOpt log_norm = lin_norm_info().YAxis(PlotOptTypes::YAxisType::log)
    .Title(PlotOptTypes::TitleType::info).LogMinimum(.2);
  PlotOpt log_norm_data = lin_norm_info().YAxis(PlotOptTypes::YAxisType::log)
    .Title(PlotOptTypes::TitleType::info)
    .LogMinimum(.2)
    .Bottom(PlotOptTypes::BottomType::ratio);
  PlotOpt lin_norm = lin_norm_info()
    .YAxis(PlotOptTypes::YAxisType::linear)
    .Title(PlotOptTypes::TitleType::info);
  PlotOpt lin_norm_data = lin_norm_info().YAxis(PlotOptTypes::YAxisType::linear)
    .Title(PlotOptTypes::TitleType::info)
    .Bottom(PlotOptTypes::BottomType::ratio);
  PlotOpt lin_norm_nooverflow = lin_norm_info()
    .YAxis(PlotOptTypes::YAxisType::linear)
    .Title(PlotOptTypes::TitleType::info)
    .Overflow(PlotOptTypes::OverflowType::none);
  PlotOpt lin_norm_nooverflow_data = lin_norm_info()
    .YAxis(PlotOptTypes::YAxisType::linear)
    .Title(PlotOptTypes::TitleType::info)
    .Overflow(PlotOptTypes::OverflowType::none)
    .Bottom(PlotOptTypes::BottomType::ratio);
  PlotOpt lin_shapes = lin_norm().Stack(PlotOptTypes::StackType::shapes).Bottom(PlotOptTypes::BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(PlotOptTypes::TitleType::info).Bottom(PlotOptTypes::BottomType::off);
  PlotOpt lin_shapes_no_overflow = lin_norm().Stack(PlotOptTypes::StackType::shapes).Bottom(PlotOptTypes::BottomType::off).Overflow(PlotOptTypes::OverflowType::none);
  PlotOpt log_shapes = lin_norm().YAxis(PlotOptTypes::YAxisType::log).Stack(PlotOptTypes::StackType::shapes).Bottom(PlotOptTypes::BottomType::ratio);
  PlotOpt log_shapes_info = lin_shapes().YAxis(PlotOptTypes::YAxisType::log).Title(PlotOptTypes::TitleType::info).Bottom(PlotOptTypes::BottomType::off);
  PlotOpt lin_lumi_shapes = lin_norm().Stack(PlotOptTypes::StackType::lumi_shapes).Title(PlotOptTypes::TitleType::info).Bottom(PlotOptTypes::BottomType::off);

  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  vector<PlotOpt> plt_lin = {lin_norm};
  vector<PlotOpt> plt_log = {log_norm};
  vector<PlotOpt> plt_shapes = {lin_shapes};
  vector<PlotOpt> plt_shapes_no_overflow = {lin_shapes_no_overflow};
  vector<PlotOpt> plt_log_shapes = {log_shapes};
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  vector<PlotOpt> plt_log_shapes_info = {log_shapes_info};
  vector<PlotOpt> plt_lin_lumi_shapes = {lin_lumi_shapes};
  if (unblind) plt_lin = {lin_norm_data};
  if (unblind) plt_log = {log_norm_data};

  // Set options
  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  //string mc_skim_folder = "mc/merged_higmc_higloose/";
  string mc_skim_folder = "mc/skim_met150/"; //needs to be loose enough to make n-1 plots, excepting MET
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  //string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt";
  //string data_skim_folder = "data/merged_higmc_higloose/";
  string data_skim_folder = "data/skim_met150/";
  string ttbar_data_skim_folder = "data/merged_higmc_higlep1T/";
  string zll_data_skim_folder = "data/merged_higmc_higlep2T/";
  string qcd_data_skim_folder = "data/merged_higmc_higqcd/";

  string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  //string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_preselect/";
  string sig_skim_folder = "SMS-TChiHH_2D/skim_met150/";
  string foldersig = mc_base_folder+year_string+"/SMS-TChiHH_2D/unskimmed/";

  //years = {2016, 2017, 2018};
  //years = {2016};
  set<int> years;
  HigUtilities::parseYears(year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(year_string);

  NamedFunc weight_notrgeff = "weight"*w_years*Functions::w_pileup;

  // Set MC 
  map<string, set<string>> mctags; 
  // Set base tags
  mctags["tt"]     = set<string>({"*TTJets_SingleLept*","TTJets_DiLept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                 "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = set<string>({"*_ST_*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root"});
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
                               "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                               "*_QCD_HT2000toInf_*",
                               "*_WH*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  vector<shared_ptr<Process> > search_signal_procs;
  //If you modify mchi_sigm, modify mchi_sigm_int in the NamedFunc below
  vector<string> mchi_sigm = {"175","500","900"}; 
  vector<string> mlsp_sigm = {"0",  "0",  "0"}; 
  vector<string> mchi_sigm2d = {"250","350","450"}; 
  vector<string> mlsp_sigm2d = {"50","200","100"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue}; // need sigm.size() <= sig_colors.size()
  vector<int> sig_colors2d = {kOrange, kYellow, kCyan};
  for (unsigned isig(0); isig<mchi_sigm.size(); isig++) {
    //search_signal_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
    //  Process::Type::signal, sig_colors[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "stitch"));
    //  std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root" << std::endl;
    search_signal_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}), "stitch"));
  }
  for (unsigned isig(0); isig<mchi_sigm2d.size(); isig++) {
    //search_signal_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
    //  Process::Type::signal, sig_colors2d[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "stitch"));
    //  std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root" << std::endl;
    search_signal_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, sig_colors2d[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}), "stitch"));
  }

  vector<shared_ptr<Process> > search_procs;
  vector<shared_ptr<Process> > cutflow_procs;
  // Set mc processes
  //search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  search_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  search_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  for (unsigned isig(0); isig<mchi_sigm.size(); isig++) {
    //search_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
    //  Process::Type::signal, sig_colors[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "stitch"));
    cutflow_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}), "stitch"));
    search_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
      Process::Type::signal, sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}), "stitch"));
    //std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root" << std::endl;
  }
  for (unsigned isig(1); isig<mchi_sigm2d.size(); isig++) {
    //search_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
    //  Process::Type::signal, sig_colors2d[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "stitch"));
    cutflow_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, sig_colors2d[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}), "stitch"));
    search_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
      Process::Type::signal, sig_colors2d[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}), "stitch"));
    //std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root" << std::endl;
  }

  if (unblind) {
    search_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                    attach_folder(data_base_folder, years, data_skim_folder, {"*.root"}),"stitch"));
  }

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  NamedFunc base_filters = Higfuncs::final_pass_filters;

  NamedFunc weight = Higfuncs::final_weight;

  NamedFunc sr_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2) && !Higfuncs::jetid_low_dphi_met &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc sr_baseline_nob = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && !Higfuncs::jetid_low_dphi_met &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc sr_baseline_nonum = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150" && !Higfuncs::jetid_low_dphi_met &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc sr_baseline_noqcd = Higfuncs::pass_filters && "met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2) &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc sr_baseline_nohig = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2) && !Higfuncs::jetid_low_dphi_met;

  NamedFunc sr_baseline_onlynum = Higfuncs::pass_filters && "met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2);

  //number loose b-jets with jet-id jet veto
  const NamedFunc jetid_nbl("jetid_nbl", [](const Baby &b) -> NamedFunc::ScalarType{
    unsigned int nbl = 0;
    if (b.pass_jets()) {
      return b.nbl();
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) continue;
      if (b.jet_deepcsv()->at(jet_idx) > 0.2217 && abs(b.SampleType()) == 2016) nbl += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.1522 && abs(b.SampleType()) == 2017) nbl += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.0494 && abs(b.SampleType()) == 2018) nbl += 1;
    }
    if (nbl > 4) nbl = 4;
    return nbl;
  });

  //number medium b-jets with jet-id jet veto
  const NamedFunc jetid_nbm("jetid_nbm", [](const Baby &b) -> NamedFunc::ScalarType{
    unsigned int nbm = 0;
    if (b.pass_jets()) {
      return b.nbm();
    }
    for (unsigned int jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
      if (!b.jet_isgood()->at(jet_idx)) continue;
      if (b.jet_id()->at(jet_idx) == false) continue;
      if (b.jet_deepcsv()->at(jet_idx) > 0.6321 && abs(b.SampleType()) == 2016) nbm += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.4941 && abs(b.SampleType()) == 2017) nbm += 1;
      if (b.jet_deepcsv()->at(jet_idx) > 0.4184 && abs(b.SampleType()) == 2018) nbm += 1;
    }
    if (nbm > 4) nbm = 4;
    return nbm;
  });

  //give the right weight for CN and N1N2 higgsino models
  const NamedFunc mixed_model_weight("mixed_model_weight",[](const Baby &b) -> NamedFunc::ScalarType{
    set<int> mchi_sigm_int = {175, 500, 900}; 
    if (b.mprod() == -999 || mchi_sigm_int.count(b.mprod()) > 0) {
      //not signal
      return Higfuncs::final_weight.GetScalar(b);
    }
    double xsec1d, xsec2d, xsec1d_unc, xsec2d_unc;
    xsec::higgsinoCrossSection(b.mprod(),xsec1d,xsec1d_unc);
    xsec::higgsino2DCrossSection(b.mprod(),xsec2d,xsec2d_unc);
    return Higfuncs::final_weight.GetScalar(b)/xsec1d*xsec2d;
  });

  // resolved cuts
  NamedFunc search_resolved = 
                         "met/mht<2 && met/met_calo<2&&weight<1.5&&"
                         "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc ttbar_resolved = 
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==1&&lep_pt[0]>30&&mt<=100&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3)||(nbt>=2&&nbm>=3&&nbl>=4))";
  NamedFunc zll_resolved =
                         "met/met_calo<5&&weight<1.5&&"
                         "nlep==2&&njet>=4&&njet<=5&&met<50&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";
  NamedFunc qcd_resolved =
                         "met/mht<2 && met/met_calo<2&&"
                         "low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
                         "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
                         "(nbm==0||nbm==1||nbm==2||nbm>=3)";

  PlotMaker pm;

  // 2b3b4b 
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved, 0, 0, weight)}), search_procs, true, true, true);

  pm.Push<Hist1D>(Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
    base_filters&&search_resolved,
    search_signal_procs, plt_log_shapes_info).Weight(weight_notrgeff).Tag("FixName:selection__search_met_signal").LuminosityTag(total_luminosity_string);
  pm.Push<Hist1D>(Axis(14, 150, 850., "met", "p_{T}^{miss} [GeV]", {200., 300., 400.}),
    base_filters&&search_resolved,
    search_procs, plt_log).Weight(weight).Tag("FixName:selection__search_met").LuminosityTag(total_luminosity_string);

  if (do_kinematic_plots) {
    //kinematic variables nb type and shape plots
    pm.Push<Hist1D>(Axis(5, 0.5, 5.5, jetid_nb, "b-tag Category (TTML)", {}),
      sr_baseline,
      search_procs, plt_log).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signalbackground_nb_ttml_log_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, 0.5, 5.5, jetid_nbm, "b-tag Category (MMMM)", {}),
      sr_baseline_nob && (jetid_nbm >= 2) && (jetid_nbm <= 4),
      search_procs, plt_log).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signalbackground_nb_mmmm_log_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(5, 0.5, 5.5, jetid_nbl, "b-tag Category (TTLL)", {}),
      sr_baseline && (jetid_nbl <= 4),
      search_procs, plt_log).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signalbackground_nb_ttll_log_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(25, 0, 250, jetid_hig_cand_am, "#LT m_{bb} #GT", {100,140}),
      Higfuncs::pass_filters && (jetid_njet >= 4) && (jetid_nb == 4),
      search_signal_procs, plt_shapes_info).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signal_am_shapes_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(30, 0, 150, jetid_hig_cand_dm, "#Delta m_{HH}", {40}),
      Higfuncs::pass_filters && (jetid_njet >= 4) && (jetid_nb == 4),
      search_signal_procs, plt_shapes_info).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signal_dm_shapes_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(40, 0, 4, jetid_hig_cand_drmax, "#Delta R_{max}", {2.2}),
      Higfuncs::pass_filters && (jetid_njet >= 4),
      search_signal_procs, plt_shapes_info).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signal_drmax_shapes_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 4, jetid_hig_cand_drmax, "#Delta R_{max}", {2.2}),
      sr_baseline_nohig && (jetid_nb == 4) && (jetid_hig_cand_am>=100 && jetid_hig_cand_am<=140),
      search_procs, plt_lin).Weight(mixed_model_weight).Tag("FixName:kinematicvars__signalbackground_drmax_lin_"+year_string).LuminosityTag(total_luminosity_string);
  }

  if (do_nminusone_plots) {
    //selection section plots
    for (int plot_type_idx = 0; plot_type_idx < 3; plot_type_idx++) {
      vector<PlotOpt> plt_type = plt_log;
      vector<PlotOpt> plt_type_signal = plt_log;
      string plt_type_string = "log";
      if (plot_type_idx == 1) {
        plt_type = plt_shapes_info;
        plt_type_signal = plt_shapes_no_overflow;
        plt_type_string = "shapes";
      }
      else if (plot_type_idx == 2) {
        plt_type = plt_lin;
        plt_type_signal = plt_lin;
        plt_type_string = "lin";
      }

      if (plot_type_idx == 1) {
        //only make shape plots right now
        //signal only plots
        pm.Push<Hist1D>(Axis(40, 150, 800, "met", "MET [GeV]", {200, 300, 450}),
          sr_baseline,
          search_signal_procs, plt_type_signal).Weight(mixed_model_weight).Tag("FixName:selection__signalcomp_met_"+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
        pm.Push<Hist1D>(Axis(20, 0, 2.2, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {1.1}),
          sr_baseline_nohig && (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200),
          search_signal_procs, plt_type_signal).Weight(mixed_model_weight).Tag("FixName:selection__signalcomp_higcanddrmax_"+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
      }

      for (int require_nb4 = 0; require_nb4 < 2; require_nb4++) {
        NamedFunc nb_cut = "1";
        string nb_desc = "";
        if (require_nb4 == 0 && plot_type_idx == 2) {
            //don't bother plotting linear plots w/o nb=4; signal too small to be visible
            continue;
        }
        if (require_nb4 == 1) {
            nb_cut = (jetid_nb==4);
            nb_desc = "nb4_";
        }
        
        // n-1 plots
        if (require_nb4==0 && plot_type_idx == 0) {
          //don't plot nb if we are requiring nb=4, right now only make log plot
          pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::jetid_nb, "N_{b}", {1.5}),
            sr_baseline_nob,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_nb_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
        }
        if (plot_type_idx == 2 && require_nb4==1) {
          //right now only make linear nb=4 plots

          pm.Push<Hist1D>(Axis(40, 150, 800, "met", "MET [GeV]", {200, 300, 450}),
            sr_baseline && nb_cut,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_met_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(8, 3.5, 11.5, Higfuncs::jetid_njet, "N_{jets}", {3.5, 5.5}),
            sr_baseline_nonum && nb_cut && "nvlep==0&&ntk==0" && (Higfuncs::jetid_nb>=2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_njet_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "nvlep", "N_{vlep}", {0.5}),
            sr_baseline_nonum && nb_cut && "ntk==0" && (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_nvlep_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(4, -0.5, 3.5, "ntk", "N_{isoTk}", {0.5}),
            sr_baseline_nonum && nb_cut && "nvlep==0" && (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_nisotk_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[0]", "#Delta#Phi_{1}", {0.5}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi1_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[1]", "#Delta#Phi_{2}", {0.5}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi2_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[2]", "#Delta#Phi_{3}", {0.3}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi3_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
          //pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[3]", "#Delta#Phi_{4}", {0.3}),
          //  HigUtilities::pass_2016 && nb_cut &&
          //  "nvlep==0 && 4 <= njet && njet <= 5 && nbt>=2 && met>150 && ntk==0 && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]<=200",
          //  search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminusphi_dphi4_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[0]", "#Delta#Phi_{1}", {0.5}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[1]>0.5 && jet_met_dphi[2]>0.3 && jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi1_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[1]", "#Delta#Phi_{2}", {0.5}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[0]>0.5 && jet_met_dphi[2]>0.3 && jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi2_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[2]", "#Delta#Phi_{3}", {0.3}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[0]>0.5 && jet_met_dphi[1]>0.5 && jet_met_dphi[3]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi3_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 3.2, "jet_met_dphi[3]", "#Delta#Phi_{4}", {0.3}),
            sr_baseline_noqcd && nb_cut && "pass_jets && jet_met_dphi[0]>0.5 && jet_met_dphi[1]>0.5 && jet_met_dphi[2]>0.3 && (met/met_calo)<2 && (met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_dphi4_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 8, "(met/mht)", "p^{miss}_{T}/H^{miss}_{T}", {2}),
            sr_baseline_noqcd && nb_cut && !Higfuncs::jetid_low_dphi_met && "(met/met_calo)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_metmht_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 8, "(met/met_calo)", "p^{miss}_{T}/p^{miss}_{Tcalo}", {2}),
            sr_baseline_noqcd && nb_cut && !Higfuncs::jetid_low_dphi_met && "(met/mht)<2",
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_metmetcalo_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 120, Higfuncs::jetid_hig_cand_dm, "#Delta m_{HH} [GeV]", {40}),
            sr_baseline_nohig && nb_cut && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_higcanddm_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 4.0, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {2.2}),
            sr_baseline_nohig && nb_cut && (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200),
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_higcanddrmax_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          pm.Push<Hist1D>(Axis(40, 0, 200, Higfuncs::jetid_hig_cand_am, "#LT m_{bb} #GT [GeV]", {100,140}),
            sr_baseline && nb_cut,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__nminus1_higcandam_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);

          //non N-1 variables
          pm.Push<Hist1D>(Axis(40, 0, 1200, "ht", "HT [GeV]", {}),
            sr_baseline && nb_cut,
            search_procs, plt_type).Weight(mixed_model_weight).Tag("FixName:selection__baseline_ht_"+nb_desc+plt_type_string+"_"+year_string).LuminosityTag(total_luminosity_string);
        }
      }
    }
  }

  if (do_cutflow) {
    pm.Push<Table>("sr_cutflow_"+year_string, vector<TableRow>{
    TableRow("Filters, $p_\\text{t}^\\text{miss}>150 \\text{ GeV}$", 
      Higfuncs::pass_filters && "met>150",0,0,mixed_model_weight),
    TableRow("$N_\\text{vl}=N_\\text{tk}=0$", 
      Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0",0,0,mixed_model_weight),
    TableRow("$4 \\leq N_\\text{jet} \\leq 5$", 
      Higfuncs::pass_filters && "met>150 && nvlep==0 && ntk==0" && (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5),0,0,mixed_model_weight),
    TableRow("$N_\\text{b}\\geq 2$", 
      sr_baseline_onlynum,0,0,mixed_model_weight),
    TableRow("Fake MET Cuts",        
      sr_baseline_nohig,0,0,mixed_model_weight),
    TableRow("$\\Delta m<40 \\text{ GeV}$, $\\langle m_\\text{bb} \\rangle<200 \\text{ GeV}$",        
      sr_baseline_nohig && (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200),0,0,mixed_model_weight),
    TableRow("$\\Delta R_\\text{max}<2.2$",        
      sr_baseline,0,0,mixed_model_weight),
    TableRow("$100 \\leq \\langle m\\rangle < 140 \\text{ GeV}$",        
      sr_baseline && (Higfuncs::jetid_hig_cand_am>=100) && (Higfuncs::jetid_hig_cand_am<140),0,0,mixed_model_weight),
    TableRow("$N_\\text{b}\\geq 3$",        
      sr_baseline && (Higfuncs::jetid_hig_cand_am>=100) && (Higfuncs::jetid_hig_cand_am<140) && (Higfuncs::jetid_nb>=3),0,0,mixed_model_weight),
    TableRow("$N_\\text{b}=4$",        
      sr_baseline && (Higfuncs::jetid_hig_cand_am>=100) && (Higfuncs::jetid_hig_cand_am<140) && (Higfuncs::jetid_nb==4),0,0,mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}>200 \\text{ GeV}$",        
      sr_baseline && (Higfuncs::jetid_hig_cand_am>=100) && (Higfuncs::jetid_hig_cand_am<140) && (Higfuncs::jetid_nb==4) && "met>200",0,0,mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}>300 \\text{ GeV}$",        
      sr_baseline && (Higfuncs::jetid_hig_cand_am>=100) && (Higfuncs::jetid_hig_cand_am<140) && (Higfuncs::jetid_nb==4) && "met>300",0,0,mixed_model_weight),
    TableRow("$p_\\text{T}^\\text{miss}>400 \\text{ GeV}$",        
      sr_baseline && (Higfuncs::jetid_hig_cand_am>=100) && (Higfuncs::jetid_hig_cand_am<140) && (Higfuncs::jetid_nb==4) && "met>400",0,0,mixed_model_weight),
    //TableRow("$\\Delta R_\\text{max}>1.1$",        
    //  "nvlep==0 && njet>=4 && njet<=5 && nbt>=2 && met>450 && ntk==0 && !low_dphi_met && (met/met_calo)<2 && (met/mht)<2 && hig_cand_dm[0]<=40 && hig_cand_drmax[0]>=1.1 && hig_cand_drmax[0]<=2.2 && hig_cand_am[0]>100 && hig_cand_am[0]<=140 && nbm>=3 && nbl>=4",0,0,mixed_model_weight),
    },cutflow_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string);
  }

  //pm.Push<Hist1D>(Axis(10,0,100,"hig_cand_dm[0]", "#Deltam [GeV]", {40.}),
  //  base_filters &&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_am[0]<200&&"
  //  "((nbt>=2&&nbm>=3&&nbl>=4))",
  //  search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_hig_cand_dm").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(10, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}),
  //  base_filters &&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_drmax[0]<2.2&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm>=3&&nbl>=4))",
  //  search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_hig_cand_am").LuminosityTag(total_luminosity_string);
  //pm.Push<Hist1D>(Axis(20,0,4,"hig_cand_drmax[0]", "#DeltaR_{max}", {1.1, 2.2}),
  //  base_filters &&
  //  "ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&"
  //  "hig_cand_am[0]<200&&hig_cand_dm[0]<40&&"
  //  "((nbt>=2&&nbm>=3&&nbl>=4))",
  //  search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_hig_cand_drmax").LuminosityTag(total_luminosity_string);

  //generate pie-charts

  //// 2b (met: 150, 200, 300, 400) low drmax
  pm.Push<Table>("FixName:selection__search_pies__2b_met150_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b_met200_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b_met300_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b_met400_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  // 3b (met: 150, 200, 300, 400) low drmax
  pm.Push<Table>("FixName:selection__search_pies__3b_met150_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__3b_met200_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__3b_met300_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__3b_met400_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  // 4b (met: 150, 200, 300, 400) low drmax
  pm.Push<Table>("FixName:selection__search_pies__4b_met150_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__4b_met200_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__4b_met300_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__4b_met400_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  // 2b (met: 150, 200, 300, 400) high drmax
  pm.Push<Table>("FixName:selection__search_pies__2b_met150_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b_met200_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b_met300_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b_met400_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt==2&&nbm==2)&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  // 3b (met: 150, 200, 300, 400) high drmax
  pm.Push<Table>("FixName:selection__search_pies__3b_met150_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__3b_met200_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__3b_met300_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__3b_met400_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm==3&&nbl==3)&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  // 4b (met: 150, 200, 300, 400) high drmax
  pm.Push<Table>("FixName:selection__search_pies__4b_met150_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__4b_met200_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__4b_met300_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__4b_met400_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"(nbt>=2&&nbm>=3&&nbl>=4)&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);

  // 2b3b4b (met: 150, 200, 300, 400) low drmax
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met150_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>150 &&met<=200 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met200_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>200 &&met<=250 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met300_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>300 &&met<=400 &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met400_lowdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>400            &&hig_cand_drmax[0]<=1.1", 0, 0, weight)}), search_procs, true, true, true);
  // 2b3b4b (met: 150, 200, 300, 400) high drmax
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met150_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>150 &&met<=200 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met200_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>200 &&met<=250 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met300_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>300 &&met<=400 &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);
  pm.Push<Table>("FixName:selection__search_pies__2b3b4b_met400_highdrmax_"+year_string  , vector<TableRow> ({TableRow("", base_filters&&search_resolved&&"nbt>=2&&met>400            &&hig_cand_drmax[0]>1.1", 0, 0, weight)}), search_procs, true, true, true);

  if (do_srambb_plots) {
    //SR <m> plots
    // 3b4b [<m>] (met: 150, 200, 300, 400) low drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met150_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met200_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met300_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]<=1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met400_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    // 3b4b [<m>] (met: 150, 200, 300, 400) high drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met150_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met200_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met300_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt>=2&&nbm>=3&&hig_cand_drmax[0]>1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_3b4b_met400_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);

    // 2b [<m>] (met: 150, 200, 300, 400) low drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met150_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met200_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met300_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]<=1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met400_lowdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    // 2b [<m>] (met: 150, 200, 300, 400) high drmax
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>150 && met<=200", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met150_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>200 && met<=300", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met200_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>300 && met<=400", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met300_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);
    pm.Push<Hist1D>(Axis(20, 0, 200, "hig_cand_am[0]", "<m_{bb}> [GeV]", {100, 140}), base_filters&&search_resolved&&"nbt==2&&nbm==2&&hig_cand_drmax[0]>1.1 && met>400            ", search_procs, plt_lin).Weight(weight).Tag("FixName:selection__search_amjj_2b_met400_highdrmax_"+year_string).LuminosityTag(total_luminosity_string);

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
      {"year", required_argument, 0, 0},
      {"unblind", no_argument, 0, 0},
      {"nokinematic", no_argument, 0, 0},
      {"nonminusone", no_argument, 0, 0},
      {"nocutflow", no_argument, 0, 0},
      {"nopiecharts", no_argument, 0, 0},
      {"nosram", no_argument, 0, 0},
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
      if(optname == "year"){
        year_string = optarg;
      } else if (optname == "unblind") {
        unblind = true;
      } else if (optname == "nokinematic") {
        do_kinematic_plots = false;
      } else if (optname == "nonminusone") {
        do_nminusone_plots = false;
      } else if (optname == "nocutflow") {
        do_cutflow = false;
      } else if (optname == "nopiecharts") {
        do_piecharts = false;
      } else if (optname == "nosram") {
        do_srambb_plots = false;
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
