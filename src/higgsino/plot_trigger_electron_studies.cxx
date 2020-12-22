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
#include "core/efficiency_plot.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"

struct ArgStruct {
  bool single_thread; // use this argument to run PlotMaker on one thread
  std::string plot_type; //"mc" or "data" 
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
  std::vector<PlotOpt> plt_lin = {lin_norm_data};
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

  std::string mc_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  std::string mc_skim_folder = "mc/skim_met150/";
  std::string sr_mc_skim_folder = "mc/merged_higmc_higloose/";
  std::string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  std::string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  std::string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  std::string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  std::string data_skim_folder = "data/skim_met150/";
  std::string sr_data_skim_folder = "data/merged_higdata_higloose/";
  std::string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  std::string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  std::string qcd_data_skim_folder = "data/merged_higdata_higqcd/";

  std::string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  std::string search_sig_skim_folder = "SMS-TChiHH_2D/skim_met150/";
  std::string sr_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";
  std::string ttbar_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep1T/";
  std::string zll_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higlep2T/";
  std::string qcd_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higqcd/";

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
  mctags["ttbar"]     = std::set<std::string>({"*TTJets_*Lept*"});
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
  mctags["noqcd"] = std::set<std::string>({"*TTJets_SingleLept*",
                               "*TTJets_DiLept*",
                               "*_TTZ*.root", "*_TTW*.root",
                               "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root", "*_ST_*.root",
                               "*_WJetsToLNu*.root", "*_ZJet*.root",
                               "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                               "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root", "*DYJetsToLL*.root"
  });

  // set qcd procs
  //std::vector<std::shared_ptr<Process> > qcd_procs;
  ////qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  ////                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ////qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  ////                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["zjets"]),"stitch"));
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["wjets"]),"stitch"));
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["single_t"]),"stitch"));
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["qcd"]),"stitch")); 
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["other"]),"stitch"));
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("Data-t#bar(t)+X", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
  //                Higfuncs::met_trigger));

  //std::vector<std::shared_ptr<Process> > qcd_procs_ttx;
  //qcd_procs_ttx.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch"));
  //qcd_procs_ttx.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
  //                Higfuncs::met_trigger));

  //// set ttbar procs
  //std::vector<std::shared_ptr<Process> > ttbar_procs;
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["zjets"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder,mctags["wjets"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["single_t"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["qcd"]),"stitch")); 
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["other"]),"stitch"));
  //ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
	//	              (Higfuncs::el_trigger || Higfuncs::mu_trigger || Higfuncs::met_trigger)));
  //
  //// set zll procs
  //std::vector<std::shared_ptr<Process> > zll_procs;
  ////zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  ////                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ////zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  ////                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["zjets"]),"stitch"));
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["wjets"]),"stitch"));
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["single_t"]),"stitch"));
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["qcd"]),"stitch")); 
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["other"]),"stitch"));
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("Data-t#bar(t)+X", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
  //                (Higfuncs::el_trigger || Higfuncs::mu_trigger)));

  //std::vector<std::shared_ptr<Process> > zll_procs_ttx;
  //zll_procs_ttx.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch"));
  //zll_procs_ttx.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
  //                attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
  //                (Higfuncs::el_trigger || Higfuncs::mu_trigger)));

  // set sr procs
  std::vector<std::shared_ptr<Process> > sr_procs;
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}", Process::Type::background, colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, {"*TTJets_SingleLeptFromT_Tune*"}), "stitch"));
  //sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", 
  //      Process::Type::background,colors("tt_1l"),{"/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/2016/mc/unskimmed/pico_TTJets_DiLept*.root"},"stitch"));

  std::vector<std::shared_ptr<Process> > mc_procs;
  mc_procs.push_back(Process::MakeShared<Baby_pico>("Real p_{T}^{miss} MC (N_{l}=0)", Process::Type::signal, kViolet,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["noqcd"]), "stitch&&nlep==0"));
  mc_procs.push_back(Process::MakeShared<Baby_pico>("Real p_{T}^{miss} MC (N_{e}=1)", Process::Type::signal, kBlue,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["noqcd"]), "stitch&&nel==1&&nlep==1"));
  mc_procs.push_back(Process::MakeShared<Baby_pico>("Real p_{T}^{miss} MC (N_{#mu}=1)", Process::Type::signal, kRed,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["noqcd"]), "stitch&&nmu==1&&nlep==1"));

  std::vector<std::shared_ptr<Process> > mc_all_procs;
  mc_all_procs.push_back(Process::MakeShared<Baby_pico>("Real p_{T}^{miss} MC (N_{l}=0)", Process::Type::signal, kViolet,
                  attach_folder(mc_base_folder, years, sr_mc_skim_folder, mctags["all"]), "stitch&&nlep==0"));
  mc_all_procs.push_back(Process::MakeShared<Baby_pico>("Real p_{T}^{miss} MC (N_{e}=1)", Process::Type::signal, kBlue,
                  attach_folder(mc_base_folder, years, sr_mc_skim_folder, mctags["all"]), "stitch&&nel==1&&nlep==1"));
  mc_all_procs.push_back(Process::MakeShared<Baby_pico>("Real p_{T}^{miss} MC (N_{#mu}=1)", Process::Type::signal, kRed,
                  attach_folder(mc_base_folder, years, sr_mc_skim_folder, mctags["all"]), "stitch&&nmu==1&&nlep==1"));

  //mc_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
  //                attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  ////mc_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
  ////                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  //mc_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
  //                attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

  std::vector<std::shared_ptr<Process> > data_procs;
  data_procs.push_back(Process::MakeShared<Baby_pico>("Data 1e", Process::Type::data, kBlue,
                  attach_folder(data_base_folder, years, data_skim_folder, {"*.root"}),
                  (Higfuncs::el_trigger && "nel==1")));
  data_procs.push_back(Process::MakeShared<Baby_pico>("Data 1#mu", Process::Type::data, kRed,
                  attach_folder(data_base_folder, years, data_skim_folder, {"*.root"}),
                  (Higfuncs::mu_trigger && "nmu==1")));

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------


  const NamedFunc met_trigger_onlymu("met_trigger_onlymu", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_met_trigger = b.HLT_PFMET90_PFMHT90_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight()||b.HLT_PFMET110_PFMHT110_IDTight()||b.HLT_PFMET120_PFMHT120_IDTight()||b.HLT_PFMET130_PFMHT130_IDTight()||b.HLT_PFMET140_PFMHT140_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight_PFHT60()||b.HLT_PFMET110_PFMHT110_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60()||b.HLT_PFMET130_PFMHT130_IDTight_PFHT60()||b.HLT_PFMET140_PFMHT140_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_HFCleaned()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned();
    return r_met_trigger;
  });

  NamedFunc sr_baseline = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150&&nvlep==0&&ntk==0" &&
    (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_njet<=5) && (Higfuncs::jetid_nb>=2) && !Higfuncs::jetid_low_dphi_met &&
    (Higfuncs::jetid_hig_cand_dm<40) && (Higfuncs::jetid_hig_cand_am<200) && (Higfuncs::jetid_hig_cand_drmax<2.2);

  NamedFunc sr_baseline_nolepcuts = Higfuncs::pass_filters && "met/met_calo<2&&met/mht<2&&met>150" &&
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

  NamedFunc weight_notrig = "weight"*Higfuncs::w_years*Functions::w_pileup;

  //------------------------------------------------------------------------------------
  //                                    plots
  //------------------------------------------------------------------------------------

  PlotMaker pm;
  if (options.plot_type == "debug") {
    pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
      Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && 
      !Higfuncs::jetid_low_dphi_met,
      Higfuncs::met_trigger,
      sr_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_mc_0el_"+options.year_string).LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
      Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && 
      !Higfuncs::jetid_low_dphi_met && "nel==1",
      Higfuncs::met_trigger,
      sr_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_mc_1el_"+options.year_string).LuminosityTag(total_luminosity_string);
  }
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && !Higfuncs::jetid_low_dphi_met,
    Higfuncs::met_trigger,
    mc_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_mc_0el_1el_1mu_"+options.year_string).LuminosityTag(total_luminosity_string);
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=4) && !Higfuncs::jetid_low_dphi_met,
    Higfuncs::met_trigger,
    mc_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_mc_0el_1el_1mu_srish_"+options.year_string).LuminosityTag(total_luminosity_string);
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && (Higfuncs::jetid_nb >= 2) && "met/mht<2&&met/met_calo<2" && !Higfuncs::jetid_low_dphi_met,
    met_trigger_onlymu,
    mc_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_mc_0el_1el_1mu_mutrigs_"+options.year_string).LuminosityTag(total_luminosity_string);
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && !Higfuncs::jetid_low_dphi_met,
    Higfuncs::met_trigger,
    mc_all_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_mc_sr_"+options.year_string).LuminosityTag(total_luminosity_string);
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && !Higfuncs::jetid_low_dphi_met,
    met_trigger_onlymu,
    data_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_data_1el_1mu_"+options.year_string).LuminosityTag(total_luminosity_string);
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=3) && !Higfuncs::jetid_low_dphi_met,
    Higfuncs::met_trigger,
    data_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_data_1el_1mu_allmettrigs_"+options.year_string).LuminosityTag(total_luminosity_string);
  pm.Push<EfficiencyPlot>(Axis(50, 150, 450, "met", "p_{T}^{miss} [GeV]", {}),
    Higfuncs::pass_filters && (Higfuncs::jetid_njet>=4) && (Higfuncs::jetid_nb == 2) && "met/mht<2&&met/met_calo<2" && !Higfuncs::jetid_low_dphi_met,
    met_trigger_onlymu,
    data_procs).Weight(weight_notrig).Tag("FixName:triggerstudies_data_1el_1mu_srish_"+options.year_string).LuminosityTag(total_luminosity_string);

  pm.Push<Table>("FixName:triggerstudies_pies_mc_sr_met150to200_comp_"+options.year_string, std::vector<TableRow> ({TableRow("", sr_baseline_nolepcuts && "150<=met&&met<200", 0, 0, weight_notrig)}), mc_all_procs, true, true, true);
  pm.Push<Table>("FixName:triggerstudies_pies_mc_sr_met150to160_comp_"+options.year_string, std::vector<TableRow> ({TableRow("", sr_baseline_nolepcuts && "150<=met&&met<160", 0, 0, weight_notrig)}), mc_all_procs, true, true, true);

  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.MakePlots(total_luminosity);
  
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
      {"plot", required_argument, 0, 0},
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
      else if (optname == "plot") {
        options.plot_type = optarg;
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

