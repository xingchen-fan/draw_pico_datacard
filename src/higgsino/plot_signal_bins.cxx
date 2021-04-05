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
}

struct ArgStruct {
  bool single_thread; // use this argument to run PlotMaker on one thread
  std::string year_string; // argument that describes years, ex. "2016", "2016,2017,2018"
  std::string out_filename; // output filename
};

//function for parsing command line arguments
ArgStruct GetOptions(int argc, char *argv[]);

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
    .Title(PlotOptTypes::TitleType::info).LogMinimum(.2).CanvasWidth(1600);
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
  std::vector<PlotOpt> plt_log = {log_norm};
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

  std::string sig_base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  std::string sig_skim_folder = "SMS-TChiHH_2D_fastSimJmeCorrection/skim_higsys/";
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

  std::vector<std::shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
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
  procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(175,0)", Process::Type::signal, kGreen,
                  attach_folder(sig_base_folder, years, sig_skim_folder, {"*mChi-175_mLSP-0_*.root"}),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(500,0)", Process::Type::signal, kRed,
                  attach_folder(sig_base_folder, years, sig_skim_folder, {"*mChi-500_mLSP-0_*.root"}),"stitch"));
  procs.push_back(Process::MakeShared<Baby_pico>("TChiHH(900,0)", Process::Type::signal, kBlue,
                  attach_folder(sig_base_folder, years, sig_skim_folder, {"*mChi-900_mLSP-0_*.root"}),"stitch"));

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------

  NamedFunc sr_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0&&njet>=4&&njet<=5&&nbt>=2&&!low_dphi_met&&hig_cand_dm[0]<40&&hig_cand_am[0]<200&&hig_cand_drmax<2.2";

  const NamedFunc signal_bin("signal_bin",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.njet() < 4) return -1;
    if (b.nbt() < 2) return -1;
    if (b.met() < 150) return -1;
    if (b.hig_cand_drmax()->at(0) > 2.2) return -1;
    int signal_bin_ = 0;
    //bin = 0 mod 2 ~ SBD, 1 mod 2 ~ HIG
    if (b.hig_cand_am()->at(0) > 100 && b.hig_cand_am()->at(0) < 140) signal_bin_ += 1;
    //bin = 0-1 mod 6 ~ 2b, 2-3 mod 6 ~ 3b, 4-5 mod 6 ~ 4b
    if (b.nbm() == 3 && b.nbl() == 3) signal_bin_ += 2;
    else if (b.nbm() >= 3 && b.nbl() >= 4) signal_bin_ += 4;
    //bin = 0-5 mod 24 ~ MET150-200, 6-11 mod 24 ~ MET200-300, 12-17 mod 24 ~ MET300-400, 18-23 mod 24 ~ MET400-Inf
    if (b.met() >= 200 && b.met() < 300) signal_bin_ += 6;
    else if (b.met() >= 300 && b.met() < 400) signal_bin_ += 12;
    else if (b.met() >= 400) signal_bin_ += 18;
    //bin = 0-23 mod 48 ~ low drmax, 24-47 mod 48 ~ high drmax
    if (b.hig_cand_drmax()->at(0) > 1.1) signal_bin_ += 24;
    return signal_bin_;
  });

  PlotMaker pm;
  int pm_idx = 0;

  //first, generate nb, drmax, and ambb plots from which to derive variations
  pm.Push<Hist1D>(Axis(48, -0.5, 47.5, signal_bin, "Signal Bin", {}),
      sr_baseline,
      procs, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:search_signalbin").LuminosityTag(total_luminosity_string);
  pm_idx++;
  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(1.0);

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

