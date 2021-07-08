//script to plot kappas after varying background composition to derive composition uncertainties

#include <unistd.h>
#include <getopt.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

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

namespace {
  //binning constants
  const std::vector<double> nb_bins = {1.5,2.5,3.5,4.5};
  const std::vector<double> drmax_bins = {0.0,1.1,2.2};
  const std::vector<double> met_bins = {150,200,300,400,600}; //includes overflow in highest bin
}

struct ArgStruct {
  bool single_thread; // use this argument to run PlotMaker on one thread
  std::string year_string; // argument that describes years, ex. "2016", "2016,2017,2018"
  std::string out_filename; // output filename
};

//function for parsing command line arguments
ArgStruct GetOptions(int argc, char *argv[]);

//function that gets variations in analysis bin variables from plot maker plots
//this requires the plot at pm_index is the nb plot, the plot at pm_index+1 is the drmax plot, and the plot at pm_index+2 is the met plot
//returns a vector of three vectors corresponding to the variations in the nb bins, drmax bins, and met bins respectively
std::vector<std::vector<double>> get_variations_from_ratio(PlotMaker &pm, unsigned int pm_index);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);
  ArgStruct options = GetOptions(argc, argv);

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
  PlotOpt log_data = lin_norm_info().YAxis(PlotOptTypes::YAxisType::log)
    .Title(PlotOptTypes::TitleType::info)
    .LogMinimum(.2)
    .Bottom(PlotOptTypes::BottomType::ratio);
  PlotOpt log_norm_data = lin_norm_info().YAxis(PlotOptTypes::YAxisType::log)
    .Title(PlotOptTypes::TitleType::info)
    .LogMinimum(.2)
    .Stack(PlotOptTypes::StackType::data_norm)
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

  std::vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};
  std::vector<PlotOpt> plt_lin = {lin_norm_info};
  std::vector<PlotOpt> plt_lin_nooverflow = {lin_norm_nooverflow_data};
  std::vector<PlotOpt> plt_log = {log_data};
  std::vector<PlotOpt> plt_log_norm = {log_norm_data};
  std::vector<PlotOpt> plt_shapes = {lin_shapes};
  std::vector<PlotOpt> plt_shapes_info = {lin_shapes_info};
  std::vector<PlotOpt> plt_lin_mc = {lin_norm};
  std::vector<PlotOpt> plt_log_mc = {log_norm};

  //------------------------------------------------------------------------------------
  //                                 samples
  //------------------------------------------------------------------------------------

  std::string mc_base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  std::string mc_skim_folder = "mc/merged_higmc_higloose/";
  std::string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  std::string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  std::string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  std::string data_base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  std::string data_skim_folder = "data/merged_higdata_higloose/";
  std::string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  std::string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  std::string qcd_data_skim_folder = "data/merged_higdata_higqcd/";

  std::string sig_base_folder = "/net/cms24/cms24r0/pico/NanoAODv7/higgsino_klamath/";
  std::string search_sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/skim_higsys/";
  //std::string search_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";
  //std::string ttbar_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep1T/";
  //std::string zll_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep2T/";
  //std::string qcd_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higqcd/";

  std::set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  float total_luminosity = 0;
  for (auto const & year : years) {
    if (year == 2016) total_luminosity += 35.9;
    if (year == 2017) total_luminosity += 41.5;
    if (year == 2018) total_luminosity += 60;
  }
  std::string total_luminosity_string = RoundNumber(total_luminosity, 1, 1).Data();

  // Set MC 
  std::map<std::string, std::set<std::string>> mctags; 
  // Set base tags
  mctags["tt"]     = std::set<std::string>({"*TTJets_*Lept*",
                                  "*_TTZ*.root", "*_TTW*.root",
                                  "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["single_t"] = std::set<std::string>({"*_ST_*.root"});
  mctags["zjets"]   = std::set<std::string>({"*_ZJet*.root", "*DYJetsToLL*.root"});
  mctags["wjets"]   = std::set<std::string>({"*_WJetsToLNu*.root"});
  mctags["qcd"]     = std::set<std::string>({"*_QCD_HT200to300_*","*_QCD_HT300to500_*","*_QCD_HT500to700_*",
                                   "*_QCD_HT700to1000_*", "*_QCD_HT1000to1500_*","*_QCD_HT1500to2000_*",
                                   "*_QCD_HT2000toInf_*"});
  mctags["other"]   = std::set<std::string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  // Combine all tags
  mctags["all"] = std::set<std::string>({"*TTJets_SingleLept*",
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

  std::vector<std::shared_ptr<Process> > signal_procs;
  std::shared_ptr<Process> signal_proc = Process::MakeShared<Baby_pico>("Signal", Process::Type::signal, kRed,
      attach_folder(sig_base_folder, years, search_sig_skim_folder, {"*_mLSP-0_*.root"}), "stitch");
  signal_procs.push_back(signal_proc);

  std::vector<std::shared_ptr<Process> > signal_data_procs;
  signal_data_procs.push_back(Process::MakeShared<Baby_pico>("(225,0)", Process::Type::signal, kRed,
      attach_folder(sig_base_folder, years, search_sig_skim_folder, {"*mChi-225_mLSP-0_*.root"}), "stitch"));
  signal_data_procs.push_back(Process::MakeShared<Baby_pico>("(400,100)", Process::Type::signal, kOrange,
      attach_folder(sig_base_folder, years, search_sig_skim_folder, {"*mChi-400_mLSP-100_*.root"}), "stitch"));
  signal_data_procs.push_back(Process::MakeShared<Baby_pico>("(400,0)", Process::Type::signal, kGreen,
      attach_folder(sig_base_folder, years, search_sig_skim_folder, {"*mChi-400_mLSP-0_*.root"}), "stitch"));
  signal_data_procs.push_back(Process::MakeShared<Baby_pico>("(800,0)", Process::Type::signal, kBlue,
      attach_folder(sig_base_folder, years, search_sig_skim_folder, {"*mChi-800_mLSP-0_*.root"}), "stitch"));
  signal_data_procs.push_back(Process::MakeShared<Baby_pico>("data", Process::Type::data, kBlack,
      attach_folder(data_base_folder, years, data_skim_folder, {"*MET*"}), Higfuncs::met_trigger));

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------

  NamedFunc sr_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2) && !Higfuncs::jetid_low_dphi_met &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc ttbar_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&nlep==1&&mt<100" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2) && (Higfuncs::lead_signal_lepton_pt>30) &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc zll_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&nlep==2&&ll_m[0]>80&&ll_m[0]<100" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::lead_signal_lepton_pt>40) &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);
  
  NamedFunc qcd_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && Higfuncs::jetid_low_dphi_met &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  const NamedFunc syst_weight_maxmurf("syst_weight_maxmurf",[](const Baby &b) -> NamedFunc::ScalarType{
    float syst_weight_maxmurf_ = 0.0;
    for (int murf_idx = 0; murf_idx < 9; murf_idx++) {
      if (b.sys_murf()->at(murf_idx) > syst_weight_maxmurf_)
        syst_weight_maxmurf_ = b.sys_murf()->at(murf_idx);
    }
    return syst_weight_maxmurf_;
  });

  const NamedFunc syst_weight_minmurf("syst_weight_minmurf",[](const Baby &b) -> NamedFunc::ScalarType{
    float syst_weight_minmurf_ = 2.0;
    for (int murf_idx = 0; murf_idx < 9; murf_idx++) {
      if (b.sys_murf()->at(murf_idx) < syst_weight_minmurf_)
        syst_weight_minmurf_ = b.sys_murf()->at(murf_idx);
    }
    return syst_weight_minmurf_;
  });

  //namedfunc for prefire since nominal weight may be zero
  const NamedFunc syst_weight_prefire_up("syst_weight_prefire_up",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.w_prefire()) < 0.00001)
      return 1.0;
    return b.sys_prefire()->at(0)/b.w_prefire();
  });

  const NamedFunc syst_weight_prefire_down("syst_weight_prefire_down",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.w_prefire()) < 0.00001)
      return 1.0;
    return b.sys_prefire()->at(1)/b.w_prefire();
  });

  const NamedFunc final_weight_nopref("final_weight_nopref",[](const Baby &b) -> NamedFunc::ScalarType{
    if (fabs(b.w_prefire()) < 0.00001)
      return 0.0;
    return Higfuncs::final_weight.GetScalar(b)/b.w_prefire();
  });

  PlotMaker pm;
  int pm_idx = 0;

  //first, generate nb, drmax, and ambb plots from which to derive variations
  pm.Push<Hist1D>(Axis(100, 0, 2.0, Higfuncs::eff_higtrig_run2_syst_up/Higfuncs::eff_higtrig_run2, 
      "Trigger upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight/Higfuncs::eff_higtrig_run2).Tag(
      "FixName:systsig_trigup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, Higfuncs::eff_higtrig_run2_syst_down/Higfuncs::eff_higtrig_run2, 
      "Trigger downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight/Higfuncs::eff_higtrig_run2).Tag(
      "FixName:systsig_trigdown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_bchig[0]/w_bhig", 
      "b-tagging upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_btagup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_bchig[1]/w_bhig", 
     "b-tagging downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_btagdown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_udsghig[0]/w_bhig", 
      "b-mistagging upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_bmisstagup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_udsghig[1]/w_bhig", 
      "b-mistagging downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_bmisstagdown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_fs_bchig[0]/w_bhig", 
      "b-tagging fastSIM upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_fsbtagup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_fs_bchig[1]/w_bhig", 
     "b-tagging fastSIM downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_fsbtagdown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_fs_udsghig[0]/w_bhig", 
      "b-mistagging fastSIM upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_fsbmisstagup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_fs_udsghig[1]/w_bhig", 
      "b-mistagging fastSIM downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_bhig").Tag(
      "FixName:systsig_fsbmisstagdown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, syst_weight_prefire_up, 
      "Prefire upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(final_weight_nopref).Tag(
      "FixName:systsig_prefireup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, syst_weight_prefire_down, 
      "Prefire downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(final_weight_nopref).Tag(
      "FixName:systsig_prefiredown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_isr[0]/w_isr", 
      "ISR upward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_isr").Tag(
      "FixName:systsig_isrup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, "sys_isr[1]/w_isr", 
      "ISR downward variations", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight*"1.0/w_isr").Tag(
      "FixName:systsig_isrdown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  //JetID, Lumi trivial
  //can't do FastSim FullSim MET, PU, JEC, JER
  pm.Push<Hist1D>(Axis(100, 0, 2.0, syst_weight_maxmurf, 
      "", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight).Tag(
      "FixName:systsig_scaleup").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, 0, 2.0, syst_weight_minmurf, 
      "", {}),
      sr_baseline,
      signal_procs, plt_lin).Weight(Higfuncs::final_weight).Tag(
      "FixName:systsig_scaledown").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.Push<Hist1D>(Axis(100, -0.5, 99.5, "npv", 
      "N_{pv}", {}),
      sr_baseline,
      signal_data_procs, plt_shapes_info).Weight(Higfuncs::final_weight).Tag(
      "FixName:systsig_pu").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(137.4);

  for (int pm_idx_2 = 0; pm_idx_2 < (pm_idx-1); pm_idx_2++) {
    Hist1D * h1d = static_cast<Hist1D*>(pm.Figures()[pm_idx_2].get());
    Hist1D::SingleHist1D* sh1d = static_cast<Hist1D::SingleHist1D*>(h1d->GetComponent(signal_proc.get()));
    TH1D hist = sh1d->scaled_hist_;
    std::cout << "Hist " << pm_idx_2 << " mean: " << hist.GetMean() << 
        ", stddev: " << hist.GetStdDev() << std::endl;
  }

  time(&endtime); 
  std::cout << std::endl << "Took " << difftime(endtime, begtime) << " seconds\n" << std::endl;
}

ArgStruct GetOptions(int argc, char *argv[]){
  ArgStruct options;
  options.single_thread = false;
  options.year_string = "2016";
  options.out_filename = "compositions_namedfunc.cpp.txt";
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"year", required_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);
    if( opt == -1) break;

    std::string optname;
    switch(opt){
    case 's':
      options.single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if (optname == "year") {
        options.year_string = optarg;
      }
      else if (optname == "out") {
        options.out_filename = optarg;
      }
      else {
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
  return options;
}

