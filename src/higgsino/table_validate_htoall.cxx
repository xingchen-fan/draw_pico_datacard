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
  string mc_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string mc_skim_folder = "mc/merged_higmc_higloose/";
  string mc_skim_folder = "mc/skim_met150/"; //needs to be loose enough to make n-1 plots, excepting MET
  string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  string data_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv5/higgsino_humboldt";
  //string data_skim_folder = "data/merged_higmc_higloose/";
  string data_skim_folder = "data/skim_met150/";
  string ttbar_data_skim_folder = "data/merged_higmc_higlep1T/";
  string zll_data_skim_folder = "data/merged_higmc_higlep2T/";
  string qcd_data_skim_folder = "data/merged_higmc_higqcd/";

  //string sig_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string sig_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/";
  string sig_fix_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath_v3/";
  //string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_preselect/";
  //string sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/skim_met150/";
  string sig_hbb_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/unskimmed/";
  string sig_hall_skim_folder = "SMS-TChiHH_HToAll_fastSimJmeCorrection/unskimmed/";
  string sig_hsep_skim_folder = "SMS-TChiHH_HToAllSeparated_fastSimJmeCorrection/unskimmed/";
  string foldersig = mc_base_folder+year_string+"/SMS-TChiHH_2D/unskimmed/";

  //string t5hh_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  //string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_preselect/";
  //string t5hh_skim_folder = "SMS-T5qqqqZH_fastSimJmeCorrection/skim_met150/";

  //string t5hhfullsim_base_folder = string(getenv("LOCAL_PICO_DIR"))+"/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string t5hh_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/";
  //string sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_preselect/";
  string t5hh_hall_skim_folder = "SMS-T5qqqqZH_FullSimJmeVariations/unskimmed/";
  string t5hh_hsep_skim_folder = "SMS-T5qqqqZH_HToAllSeparated_FullSimJmeVariations/unskimmed/";

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
  //vector<string> mchi_sigm = {"600","800","1000"}; 
  //vector<string> mlsp_sigm = {"0",  "0",  "0"}; 
  //vector<string> mchi_sigm2d = {"250","350","450"}; 
  //vector<string> mlsp_sigm2d = {"50","200","100"}; 
  //vector<string> t5hh_filenames = {"*mGluino-1600_mChi-1550_mLSP-1_*.root","*mGluino-2000_mChi-1950_mLSP-1_*.root","*mGluino-2200_mChi-2150_mLSP-1_*.root"}; 
  //vector<string> t5hh_modelnames = {"T5HH1600","T5HH2000","T5HH2200"}; 
  //vector<string> t5hhfullsim_modelnames = {"T5HH1600(FullSIM)","T5HH2000(FullSIM)","T5HH2200(FullSIM)"}; 
  //vector<int> sig_colors = {kGreen+1, kRed, kBlue}; // need sigm.size() <= sig_colors.size()
  //vector<int> sig_colors2d = {kOrange, kYellow, kCyan};
  //vector<int> sig_colorst5hh = {kOrange, kYellow, kCyan};

  vector<string> sig_name_string = {"(400,1)","(600,1)"};
  vector<string> sig_file_string = {"mChi-400_mLSP-0","mChi-600_mLSP-0"};
  vector<int> sig_hbb_colors = {kRed,kBlue};
  vector<int> sig_hall_colors = {kRed,kBlue};
  vector<int> sig_hsep_colors = {kRed,kBlue};

  vector<string> t5hh_name_string = {"(1000,1)","(1600,1)","(2000,1)"};
  vector<string> t5hh_file_string = {"mGluino-1000_mChi-950_mLSP-1","mGluino-1600_mChi-1550_mLSP-1","mGluino-2000_mChi-1950_mLSP-1"};
  vector<int> t5hh_hbb_colors = {kRed,kBlue,kGreen};
  vector<int> t5hh_hall_colors = {kRed,kBlue,kGreen};
  vector<int> t5hh_hsep_colors = {kRed,kBlue,kGreen};

  vector<shared_ptr<Process> > cutflow_procs;
  for (unsigned isig(0); isig<sig_name_string.size(); isig++) {
    cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH"+sig_name_string[isig]+" (Pre-fix)", 
      Process::Type::signal, sig_hbb_colors[isig], attach_folder(sig_base_folder, years, sig_hbb_skim_folder, {"*TChiHH*_"+sig_file_string[isig]+"*.root"}), "stitch"));
    cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH"+sig_name_string[isig]+" (Post-fix)", 
      Process::Type::signal, sig_hall_colors[isig], attach_folder(sig_fix_base_folder, years, sig_hbb_skim_folder, {"*TChiHH*_"+sig_file_string[isig]+"*.root"}), "stitch"));
    //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH"+sig_name_string[isig]+" (H to all,separate)", 
    //  Process::Type::signal, sig_hsep_colors[isig], attach_folder(sig_base_folder, years, sig_hsep_skim_folder, {"*TChiHH*_"+sig_file_string[isig]+"*.root"}), "stitch"));
  }

  //for (unsigned isig(0); isig<t5hh_name_string.size(); isig++) {
  //  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("T5HH"+t5hh_name_string[isig]+" (H to all,combined)", 
  //    Process::Type::signal, t5hh_hall_colors[isig], attach_folder(t5hh_base_folder, years, t5hh_hall_skim_folder, {"*"+t5hh_file_string[isig]+"*.root"}), "stitch"));
  //  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("T5HH"+t5hh_name_string[isig]+" (H to all,separate)", 
  //    Process::Type::signal, t5hh_hsep_colors[isig], attach_folder(t5hh_base_folder, years, t5hh_hsep_skim_folder, {"*"+t5hh_file_string[isig]+"*.root"}), "stitch"));
  //    //Process::Type::signal, t5hh_hsep_colors[isig], attach_folder(t5hh_base_folder, years, t5hh_hsep_skim_folder, {"*HToBB_HToBB*"+t5hh_file_string[isig]+"*.root", "*HToBB_HToWW*"+t5hh_file_string[isig]+"*.root"}), "stitch"));
  //    //Process::Type::signal, t5hh_hsep_colors[isig], attach_folder(t5hh_base_folder, years, t5hh_hsep_skim_folder, {"*HToBB_HToBB*"+t5hh_file_string[isig]+"*.root", "*HToBB_HToWW*"+t5hh_file_string[isig]+"*.root", "*HToBB_HToGluGlu*"+t5hh_file_string[isig]+"*.root", "*HToBB_HToTauTau*"+t5hh_file_string[isig]+"*.root", "*HToBB_HToCC*"+t5hh_file_string[isig]+"*.root", "*HToBB_HToZZ*"+t5hh_file_string[isig]+"*.root"}), "stitch"));
  //}

  //// Set mc processes
  //////search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //////                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //////search_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //////                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  //search_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  //search_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch"));
  //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("V+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["vjets"]),"stitch"));
  //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  //for (unsigned isig(0); isig<mchi_sigm.size(); isig++) {
  //  //search_procs.push_back(Process::MakeShared<Baby_pico>("GMSB("+mchi_sigm[isig]+")", 
  //  //  Process::Type::signal, sig_colors[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}, "stitch"));
  //  cutflow_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH"+mchi_sigm[isig], 
  //    Process::Type::signal, sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}), "stitch"));
  //  //search_procs.push_back(Process::MakeShared<Baby_pico>("TChiHH"+mchi_sigm[isig], 
  //  //  Process::Type::signal, sig_colors[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root"}), "stitch"));
  //  //std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm[isig]+"_mLSP-"+mlsp_sigm[isig]+"*.root" << std::endl;
  //}
  //for (unsigned isig(1); isig<mchi_sigm2d.size(); isig++) {
  //  //search_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
  //  //  Process::Type::signal, sig_colors2d[isig], {foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}, "stitch"));
  //  //cutflow_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
  //  //  Process::Type::signal, sig_colors2d[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}), "stitch"));
  //  //search_procs.push_back(Process::MakeShared<Baby_pico>("("+mchi_sigm2d[isig]+","+mlsp_sigm2d[isig]+")", 
  //  //  Process::Type::signal, sig_colors2d[isig], attach_folder(sig_base_folder, years, sig_skim_folder, {"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root"}), "stitch"));
  //  //std::cout << foldersig+"*TChiHH_mChi-"+mchi_sigm2d[isig]+"_mLSP-"+mlsp_sigm2d[isig]+"*.root" << std::endl;
  //}
  ////for (unsigned isig(0); isig<t5hh_filenames.size(); isig++) {
  ////  cutflow_procs.push_back(Process::MakeShared<Baby_pico>(t5hh_modelnames[isig], 
  ////    Process::Type::signal, sig_colorst5hh[isig], attach_folder(t5hh_base_folder, years, t5hh_skim_folder, {t5hh_filenames[isig]}), "stitch"));
  ////}
  //for (unsigned isig(0); isig<t5hh_filenames.size(); isig++) {
  //  cutflow_procs.push_back(Process::MakeShared<Baby_pico>(t5hh_modelnames[isig], 
  //    Process::Type::signal, sig_colorst5hh[isig], attach_folder(t5hh_base_folder, years, t5hh_skim_folder, {t5hh_filenames[isig]}), "stitch"));
  //}

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------
  
  NamedFunc base_filters = Higfuncs::final_pass_filters;

  const NamedFunc anti_hhbbbb_br("anti_hhbbbb_br", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.type()/1000 == 107) return 1.0/0.56/0.56;
    return 1.0;
  });

  NamedFunc weight = Higfuncs::final_weight*anti_hhbbbb_br;

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

  const NamedFunc hhtobbbb("hhtobbbb", [](const Baby &b) -> NamedFunc::ScalarType{
    if ((b.type() / 1000 != 107) && (b.type() / 1000 != 106)) return true; //non SUSY model
    int num_bs = 0;
    for (unsigned int mc_idx(0); mc_idx < b.mc_id()->size(); mc_idx++) {
      if (abs(b.mc_id()->at(mc_idx)) == 5 && b.mc_mom()->at(mc_idx) == 25)
        num_bs++;
    }
    if (num_bs==4) return true;
    return false;
  });

  const NamedFunc hhtobbww("hhtobbww", [](const Baby &b) -> NamedFunc::ScalarType{
    if ((b.type() / 1000 != 107) && (b.type() / 1000 != 106)) return true; //non SUSY model
    int num_bs = 0;
    int num_ws = 0;
    for (unsigned int mc_idx(0); mc_idx < b.mc_id()->size(); mc_idx++) {
      if (abs(b.mc_id()->at(mc_idx)) == 5 && b.mc_mom()->at(mc_idx) == 25)
        num_bs++;
      if (abs(b.mc_id()->at(mc_idx)) == 24 && b.mc_mom()->at(mc_idx) == 25)
        num_ws++;
    }
    //if (num_bs == 4) return true;
    if (num_bs==2&&num_ws==2) return true;
    return false;
  });

  const NamedFunc hhtobbgluglu("hhtobbgluglu", [](const Baby &b) -> NamedFunc::ScalarType{
    if ((b.type() / 1000 != 107) && (b.type() / 1000 != 106)) return true; //non SUSY model
    int num_bs = 0;
    int num_glus = 0;
    for (unsigned int mc_idx(0); mc_idx < b.mc_id()->size(); mc_idx++) {
      if (abs(b.mc_id()->at(mc_idx)) == 5 && b.mc_mom()->at(mc_idx) == 25)
        num_bs++;
      if (abs(b.mc_id()->at(mc_idx)) == 21 && b.mc_mom()->at(mc_idx) == 25)
        num_glus++;
    }
    //if (num_bs == 4) return true;
    if (num_bs==2&&num_glus==2) return true;
    return false;
  });

  const NamedFunc hhtobbtautau("hhtobbtautau", [](const Baby &b) -> NamedFunc::ScalarType{
    if ((b.type() / 1000 != 107) && (b.type() / 1000 != 106)) return true; //non SUSY model
    int num_bs = 0;
    int num_taus = 0;
    for (unsigned int mc_idx(0); mc_idx < b.mc_id()->size(); mc_idx++) {
      if (abs(b.mc_id()->at(mc_idx)) == 5 && b.mc_mom()->at(mc_idx) == 25)
        num_bs++;
      if (abs(b.mc_id()->at(mc_idx)) == 15 && b.mc_mom()->at(mc_idx) == 25)
        num_taus++;
    }
    //if (num_bs == 4) return true;
    if (num_bs==2&&num_taus==2) return true;
    return false;
  });

  const NamedFunc hhtobbcc("hhtobbcc", [](const Baby &b) -> NamedFunc::ScalarType{
    if ((b.type() / 1000 != 107) && (b.type() / 1000 != 106)) return true; //non SUSY model
    int num_bs = 0;
    int num_cs = 0;
    for (unsigned int mc_idx(0); mc_idx < b.mc_id()->size(); mc_idx++) {
      if (abs(b.mc_id()->at(mc_idx)) == 5 && b.mc_mom()->at(mc_idx) == 25)
        num_bs++;
      if (abs(b.mc_id()->at(mc_idx)) == 4 && b.mc_mom()->at(mc_idx) == 25)
        num_cs++;
    }
    //if (num_bs == 4) return true;
    if (num_bs==2&&num_cs==2) return true;
    return false;
  });

  const NamedFunc hhtobbzz("hhtobbzz", [](const Baby &b) -> NamedFunc::ScalarType{
    if ((b.type() / 1000 != 107) && (b.type() / 1000 != 106)) return true; //non SUSY model
    int num_bs = 0;
    int num_zs = 0;
    for (unsigned int mc_idx(0); mc_idx < b.mc_id()->size(); mc_idx++) {
      if (abs(b.mc_id()->at(mc_idx)) == 5 && b.mc_mom()->at(mc_idx) == 25)
        num_bs++;
      if (abs(b.mc_id()->at(mc_idx)) == 23 && b.mc_mom()->at(mc_idx) == 25)
        num_zs++;
    }
    //if (num_bs == 4) return true;
    if (num_bs==2&&num_zs==2) return true;
    return false;
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

  pm.Push<Table>("validate_htoall_"+year_string, vector<TableRow>{
  //TableRow("No Cuts (w\\_lumi)", 
  //  "1",0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ bbbb (w\\_lumi)", 
  //  hhtobbbb,0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ bbWW (w\\_lumi)", 
  //  hhtobbww,0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ bbgg (w\\_lumi)", 
  //  hhtobbgluglu,0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ bb$\\tau\\tau$ (w\\_lumi)", 
  //  hhtobbtautau,0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ bbcc (w\\_lumi)", 
  //  hhtobbcc,0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ bbZZ (w\\_lumi)", 
  //  hhtobbzz,0,0,"w_lumi"*w_years_search),

  //TableRow("HH$\\to$ others (w\\_lumi)", 
  //  (!hhtobbbb)&&(!hhtobbww)&&(!hhtobbgluglu)&&(!hhtobbtautau)&&(!hhtobbcc)&&(!hhtobbzz),0,0,"w_lumi"*w_years_search),

  TableRow("H$\\to$ bb (w\\_lumi)", 
    Higfuncs::htobb,0,0,"w_lumi"*w_years_search),

  TableRow("H$\\to$ bb (lep,lumi,isr weights)", 
    Higfuncs::htobb,0,0,"w_lumi*w_lep*w_fs_lep"*w_years_search),

  TableRow("H$\\to$ bb (lep,lumi,isr,btag weights)", 
    Higfuncs::htobb,0,0,"w_lumi*w_lep*w_fs_lep*w_bhig"*w_years_search),

  TableRow("H$\\to$ bb (lep,lumi,isr,btag,pf weights)", 
    Higfuncs::htobb,0,0,"w_lumi*w_lep*w_fs_lep*w_bhig*w_prefire"*w_years_search),

  TableRow("H$\\to$ bb (all weights)", 
    Higfuncs::htobb,0,0,Higfuncs::final_weight),

  TableRow("$p_\\text{T}^\\text{miss}>150$ GeV", 
    Higfuncs::htobb && "met>150",0,0,Higfuncs::final_weight),

  TableRow("$4\\leq N_\\text{j} \\leq 5$", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5",0,0,Higfuncs::final_weight),

  TableRow("$N_\\text{l}=N_\\text{tk}=0$", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0",0,0,Higfuncs::final_weight),

  TableRow("$N_\\text{bT}\\geq 2$", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2",0,0,Higfuncs::final_weight),

  TableRow("$p_\\text{T}^\\text{miss}$ quality", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met",0,0,Higfuncs::final_weight),

  TableRow("$\\Delta m < 40$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40",0,0,Higfuncs::final_weight),

  TableRow("$\\langle m \\rangle < 200$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200",0,0,Higfuncs::final_weight),

  TableRow("$\\Delta R_\\text{max} < 2.2$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2",0,0,Higfuncs::final_weight),

  TableRow("SR 3b-High $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}<300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm==3&&nbl==3&&hig_cand_drmax[0]>1.1&&met<300",0,0,Higfuncs::final_weight),

  TableRow("SR 3b-High $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}>300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm==3&&nbl==3&&hig_cand_drmax[0]>1.1&&met>300",0,0,Higfuncs::final_weight),

  TableRow("SR 3b-Low $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}<300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm==3&&nbl==3&&hig_cand_drmax[0]<1.1&&met<300",0,0,Higfuncs::final_weight),

  TableRow("SR 3b-Low $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}>300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm==3&&nbl==3&&hig_cand_drmax[0]<1.1&&met>300",0,0,Higfuncs::final_weight),

  TableRow("SR 4b-High $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}<300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm>=3&&nbl>=4&&hig_cand_drmax[0]>1.1&&met<300",0,0,Higfuncs::final_weight),

  TableRow("SR 4b-High $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}>300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm>=3&&nbl>=4&&hig_cand_drmax[0]>1.1&&met>300",0,0,Higfuncs::final_weight),

  TableRow("SR 4b-Low $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}<300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm>=3&&nbl>=4&&hig_cand_drmax[0]<1.1&&met<300",0,0,Higfuncs::final_weight),

  TableRow("SR 4b-Low $\\Delta R_\\text{max}$-$p_\\text{T}^\\text{miss}>300$ GeV", 
    Higfuncs::htobb && "met>150&&4<=njet&&njet<=5&&nlep==0&&ntk==0&&nbt>=2&&(met/mht<2)&&(met/met_calo<2)&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax[0]<2.2&&nbm>=3&&nbl>=4&&hig_cand_drmax[0]<1.1&&met>300",0,0,Higfuncs::final_weight),

  },cutflow_procs,false,true,false,true,false,true).LuminosityTag(total_luminosity_string);

  //pm.Push<Table>("FixName:boo_baseline_pie"  , vector<TableRow> ({TableRow("", 
  //    Higfuncs::pass_filters && "njet>=1&&met>300&&mht>=200&&ntk==0&&nvlep==0&&!low_dphi_met&&ht>600" && (nfatjet>=2) && fatjet_ptcut && fatjet_masscut_loose,0,0,Higfuncs::final_weight)
  //    }), search_procs, true, true, true);

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
