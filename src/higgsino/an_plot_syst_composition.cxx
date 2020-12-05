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
  std::string mc_skim_folder = "mc/merged_higmc_higloose/";
  std::string ttbar_mc_skim_folder = "mc/merged_higmc_higlep1T/";
  std::string zll_mc_skim_folder = "mc/merged_higmc_higlep2T/";
  std::string qcd_mc_skim_folder = "mc/merged_higmc_higqcd/";

  std::string data_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  std::string data_skim_folder = "data/merged_higdata_higloose/";
  std::string ttbar_data_skim_folder = "data/merged_higdata_higlep1T/";
  std::string zll_data_skim_folder = "data/merged_higdata_higlep2T/";
  std::string qcd_data_skim_folder = "data/merged_higdata_higqcd/";

  std::string sig_base_folder = "/net/cms25/cms25r5/pico/NanoAODv7/higgsino_inyo/";
  std::string search_sig_skim_folder = "SMS-TChiHH_2D/merged_higmc_higloose/";
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

  // set qcd procs
  std::vector<std::shared_ptr<Process> > qcd_procs;
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //qcd_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["zjets"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder,mctags["wjets"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["single_t"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["qcd"]),"stitch")); 
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["other"]),"stitch"));
  qcd_procs.push_back(Process::MakeShared<Baby_pico>("Data-t#bar(t)+X", Process::Type::data, kBlack,
                  attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
                  Higfuncs::met_trigger));

  std::vector<std::shared_ptr<Process> > qcd_procs_ttx;
  qcd_procs_ttx.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, qcd_mc_skim_folder, mctags["tt"]),"stitch"));
  qcd_procs_ttx.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                  attach_folder(data_base_folder, years, qcd_data_skim_folder, {"*.root"}),
                  Higfuncs::met_trigger));

  // set ttbar procs
  std::vector<std::shared_ptr<Process> > ttbar_procs;
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, ttbar_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
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
  ttbar_procs.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                  attach_folder(data_base_folder, years, ttbar_data_skim_folder, {"*.root"}),
		              (Higfuncs::el_trigger || Higfuncs::mu_trigger || Higfuncs::met_trigger)));
  
  // set zll procs
  std::vector<std::shared_ptr<Process> > zll_procs;
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  //zll_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
  //                attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["zjets"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder,mctags["wjets"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["single_t"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["qcd"]),"stitch")); 
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["other"]),"stitch"));
  zll_procs.push_back(Process::MakeShared<Baby_pico>("Data-t#bar(t)+X", Process::Type::data, kBlack,
                  attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                  (Higfuncs::el_trigger || Higfuncs::mu_trigger)));

  std::vector<std::shared_ptr<Process> > zll_procs_ttx;
  zll_procs_ttx.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, zll_mc_skim_folder, mctags["tt"]),"stitch"));
  zll_procs_ttx.push_back(Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack,
                  attach_folder(data_base_folder, years, zll_data_skim_folder, {"*.root"}),
                  (Higfuncs::el_trigger || Higfuncs::mu_trigger)));

  // set sr procs
  std::vector<std::shared_ptr<Process> > sr_procs;
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}>0)", Process::Type::background,colors("tt_htau"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh>0"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("t#bar{t}+X (#tau_{had}=0)", Process::Type::background,colors("tt_1l"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["tt"]),"stitch&&ntrutauh==0"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Z+jets", Process::Type::background, kOrange+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["zjets"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("W+jets", Process::Type::background, kGreen+1,
                  attach_folder(mc_base_folder, years, mc_skim_folder,mctags["wjets"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Single t", Process::Type::background,colors("single_t"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["single_t"]),"stitch"));
  sr_procs.push_back(Process::MakeShared<Baby_pico>("QCD", Process::Type::background, colors("other"),
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["qcd"]),"stitch")); 
  sr_procs.push_back(Process::MakeShared<Baby_pico>("Other", Process::Type::background, kGray+2,
                  attach_folder(mc_base_folder, years, mc_skim_folder, mctags["other"]),"stitch"));

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

  //------------------------------------------------------------------------------------
  //                                   ttX subtraction plots
  //------------------------------------------------------------------------------------

  PlotMaker pm_ttx_subtraction;
  pm_ttx_subtraction.Push<Hist1D>(Axis(nb_bins, Higfuncs::jetid_nb, "N_{b}", {}),
    zll_baseline && Higfuncs::jetid_nb>=2,
    zll_procs_ttx, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttxsub_zll_nb").LuminosityTag(total_luminosity_string);
  pm_ttx_subtraction.Push<Hist1D>(Axis(drmax_bins, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {}),
    zll_baseline,
    zll_procs_ttx, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttxsub_zll_drmax").LuminosityTag(total_luminosity_string);
  pm_ttx_subtraction.Push<Hist1D>(Axis(met_bins, "ll_pt[0]", "p_{Tll} [GeV]", {}),
    zll_baseline && "ll_pt[0]>=150",
    zll_procs_ttx, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttxsub_zll_met").LuminosityTag(total_luminosity_string);
  pm_ttx_subtraction.Push<Hist1D>(Axis(nb_bins, Higfuncs::jetid_nb, "N_{b}", {}),
    qcd_baseline && Higfuncs::jetid_nb>=2,
    qcd_procs_ttx, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttxsub_qcd_nb").LuminosityTag(total_luminosity_string);
  pm_ttx_subtraction.Push<Hist1D>(Axis(drmax_bins, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {}),
    qcd_baseline,
    qcd_procs_ttx, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttxsub_qcd_drmax").LuminosityTag(total_luminosity_string);
  pm_ttx_subtraction.Push<Hist1D>(Axis(met_bins, "met", "p_{T}^{miss} [GeV]", {}),
    qcd_baseline,
    qcd_procs_ttx, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttxsub_qcd_met").LuminosityTag(total_luminosity_string);

  pm_ttx_subtraction.multithreaded_ = !options.single_thread;
  pm_ttx_subtraction.min_print_ = true;
  pm_ttx_subtraction.MakePlots(1.0);

  //------------------------------------------------------------------------------------
  //                              make data scaling NamedFuncs
  //------------------------------------------------------------------------------------
  
  std::vector<std::vector<double>> zll_ttx_subtraction = get_variations_from_ratio(pm_ttx_subtraction, 0);
  std::vector<std::vector<double>> qcd_ttx_subtraction = get_variations_from_ratio(pm_ttx_subtraction, 3);

  //return a weight that removes tt+X from data for nb plots
  const NamedFunc ttx_subtraction_weight_nb("ttx_subtraction_weight_nb",[zll_ttx_subtraction, qcd_ttx_subtraction](const Baby &b) -> NamedFunc::ScalarType {
    double weight = Higfuncs::final_weight.GetScalar(b);
    if (b.SampleType() > 0) return weight; //MC
    if (b.nlep() == 1 || (b.nlep() == 0 && !b.low_dphi_met())) return weight; //1l CR/SR
    for (unsigned int nb_bin = 0; nb_bin < nb_bins.size()-1; nb_bin++) {
      if (Higfuncs::jetid_nb.GetScalar(b) > nb_bins[nb_bin] && Higfuncs::jetid_nb.GetScalar(b) < nb_bins[nb_bin+1] ) {
        if (b.nlep() == 0) weight = weight*(1.0-1.0/qcd_ttx_subtraction[0][nb_bin]);
        if (b.nlep() == 2) weight = weight*(1.0-1.0/zll_ttx_subtraction[0][nb_bin]);
      }
    }
    return weight;
  });

  //return a weight that removes tt+X from data for drmax plots
  const NamedFunc ttx_subtraction_weight_drmax("ttx_subtraction_weight_drmax",[zll_ttx_subtraction, qcd_ttx_subtraction](const Baby &b) -> NamedFunc::ScalarType {
    double weight = Higfuncs::final_weight.GetScalar(b);
    if (b.SampleType() > 0) return weight; //MC
    if (b.nlep() == 1 || (b.nlep() == 0 && !b.low_dphi_met())) return weight; //1l CR/SR
    for (unsigned int drmax_bin = 0; drmax_bin < drmax_bins.size()-1; drmax_bin++) {
      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > drmax_bins[drmax_bin] && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < drmax_bins[drmax_bin+1] ) {
        if (b.nlep() == 0) weight = weight*(1.0-1.0/qcd_ttx_subtraction[1][drmax_bin]);
        if (b.nlep() == 2) weight = weight*(1.0-1.0/zll_ttx_subtraction[1][drmax_bin]);
      }
    }
    return weight;
  });

  //return a weight that removes tt+X from data for met plots
  const NamedFunc ttx_subtraction_weight_met("ttx_subtraction_weight_met",[zll_ttx_subtraction, qcd_ttx_subtraction](const Baby &b) -> NamedFunc::ScalarType {
    double weight = Higfuncs::final_weight.GetScalar(b);
    if (b.SampleType() > 0) return weight; //MC
    if (b.nlep() == 1 || (b.nlep() == 0 && !b.low_dphi_met())) return weight; //1l CR/SR
    for (unsigned int met_bin = 0; met_bin < met_bins.size()-1; met_bin++) {
      //set the upper bound to the edge of the bin or 9999 if top bin
      double upper_bound = met_bins[met_bin+1];
      if (met_bin == met_bins.size()-2) {
        upper_bound = 9999;
      }
      if (b.met() > met_bins[met_bin] && b.met() < upper_bound ) {
        if (b.nlep() == 0) weight = weight*(1.0-1.0/qcd_ttx_subtraction[2][met_bin]);
        if (b.nlep() == 2) weight = weight*(1.0-1.0/zll_ttx_subtraction[2][met_bin]);
      }
    }
    return weight;
  });

  //------------------------------------------------------------------------------------
  //                                   Variation Plots
  //------------------------------------------------------------------------------------

  PlotMaker pm_variations;

  //first, generate nb, drmax, and ambb plots from which to derive variations
  pm_variations.Push<Hist1D>(Axis(nb_bins, Higfuncs::jetid_nb, "N_{b}", {}),
    ttbar_baseline,
    ttbar_procs, plt_log_norm).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttbar_nb").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(drmax_bins, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {}),
    ttbar_baseline,
    ttbar_procs, plt_log_norm).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttbar_drmax").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(met_bins, "met", "p_{T}^{miss} [GeV]", {}),
    ttbar_baseline && "met>=150",
    ttbar_procs, plt_log).Weight(Higfuncs::final_weight).Tag("FixName:syst_composition_ttbar_met").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(nb_bins, Higfuncs::jetid_nb, "N_{b}", {}),
    zll_baseline && Higfuncs::jetid_nb>=2,
    zll_procs, plt_log_norm).Weight(ttx_subtraction_weight_nb).Tag("FixName:syst_composition_zll_nb").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(drmax_bins, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {}),
    zll_baseline,
    zll_procs, plt_log_norm).Weight(ttx_subtraction_weight_drmax).Tag("FixName:syst_composition_zll_drmax").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(met_bins, "ll_pt[0]", "p_{Tll} [GeV]", {}),
    zll_baseline && "ll_pt[0]>=150",
    zll_procs, plt_log).Weight(ttx_subtraction_weight_met).Tag("FixName:syst_composition_zll_met").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(nb_bins, Higfuncs::jetid_nb, "N_{b}", {}),
    qcd_baseline && Higfuncs::jetid_nb>=2,
    qcd_procs, plt_log_norm).Weight(ttx_subtraction_weight_nb).Tag("FixName:syst_composition_qcd_nb").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(drmax_bins, Higfuncs::jetid_hig_cand_drmax, "#Delta R_{max}", {}),
    qcd_baseline,
    qcd_procs, plt_log_norm).Weight(ttx_subtraction_weight_drmax).Tag("FixName:syst_composition_qcd_drmax").LuminosityTag(total_luminosity_string);
  pm_variations.Push<Hist1D>(Axis(met_bins, "met", "p_{T}^{miss} [GeV]", {}),
    qcd_baseline,
    qcd_procs, plt_log).Weight(ttx_subtraction_weight_met).Tag("FixName:syst_composition_qcd_met").LuminosityTag(total_luminosity_string);
  
  pm_variations.multithreaded_ = !options.single_thread;
  pm_variations.min_print_ = true;
  pm_variations.MakePlots(1.0);

  //------------------------------------------------------------------------------------
  //                           Get Variations From Ratio Plots
  //------------------------------------------------------------------------------------
  
  //get the histograms from the plot maker
  std::vector<std::vector<double>> ttbar_variations = get_variations_from_ratio(pm_variations, 0);
  std::vector<std::vector<double>> zll_variations = get_variations_from_ratio(pm_variations, 3);
  std::vector<std::vector<double>> qcd_variations = get_variations_from_ratio(pm_variations, 6);
  
  std::ofstream comp_output_file;
  comp_output_file.open(("src/higgsino/"+options.out_filename).c_str());
  comp_output_file << "//return a weight that can be multiplied with particular MC samples for varying composition;\n";
  comp_output_file << "const NamedFunc composition_variation_weight(\"composition_variation_weight\",[](const Baby &b) -> NamedFunc::ScalarType{\n";
  comp_output_file << "  double weight = 1.0;\n";
  comp_output_file << "  int sample = b.type()/1000;\n";
  comp_output_file << "  switch (sample) {\n";
  comp_output_file << "    case 1:\n";
  bool first_bin = true;
  for (unsigned int nb_bin = 0; nb_bin < nb_bins.size()-1; nb_bin++) {
    if (first_bin) {
      comp_output_file << "      if (Higfuncs::jetid_nb.GetScalar(b) > " << nb_bins[nb_bin] << " && Higfuncs::jetid_nb.GetScalar(b) < " << nb_bins[nb_bin+1] << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (Higfuncs::jetid_nb.GetScalar(b) > " << nb_bins[nb_bin] << " && Higfuncs::jetid_nb.GetScalar(b) < " << nb_bins[nb_bin+1] << ") \n";
    }
    comp_output_file << "        weight = weight*" << ttbar_variations[0][nb_bin] <<";\n";
  }
  first_bin = true;
  for (unsigned int drmax_bin = 0; drmax_bin < drmax_bins.size()-1; drmax_bin++) {
    if (first_bin) {
      comp_output_file << "      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > " << drmax_bins[drmax_bin] << " && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < " << drmax_bins[drmax_bin+1] << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > " << drmax_bins[drmax_bin] << " && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < " << drmax_bins[drmax_bin+1] << ") \n";
    }
    comp_output_file << "        weight = weight*" << ttbar_variations[1][drmax_bin] <<";\n";
  }
  first_bin = true;
  for (unsigned int met_bin = 0; met_bin < met_bins.size()-1; met_bin++) {
    //set the upper bound to the edge of the bin or 9999 if top bin
    double upper_bound = met_bins[met_bin+1];
    if (met_bin == met_bins.size()-2) {
      upper_bound = 9999;
    }
    if (first_bin) {
      comp_output_file << "      if (b.met() > " << met_bins[met_bin] << " && b.met() < " << upper_bound << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (b.met() > " << met_bins[met_bin] << " && b.met() < " << upper_bound << ") \n";
    }
    comp_output_file << "        weight = weight*" << ttbar_variations[2][met_bin] <<";\n";
  }
  comp_output_file << "      break;\n";
  comp_output_file << "    case 8:\n";
  first_bin = true;
  for (unsigned int nb_bin = 0; nb_bin < nb_bins.size()-1; nb_bin++) {
    if (first_bin) {
      comp_output_file << "      if (Higfuncs::jetid_nb.GetScalar(b) > " << nb_bins[nb_bin] << " && Higfuncs::jetid_nb.GetScalar(b) < " << nb_bins[nb_bin+1] << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (Higfuncs::jetid_nb.GetScalar(b) > " << nb_bins[nb_bin] << " && Higfuncs::jetid_nb.GetScalar(b) < " << nb_bins[nb_bin+1] << ") \n";
    }
    comp_output_file << "        weight = weight*" << zll_variations[0][nb_bin] <<";\n";
  }
  first_bin = true;
  for (unsigned int drmax_bin = 0; drmax_bin < drmax_bins.size()-1; drmax_bin++) {
    if (first_bin) {
      comp_output_file << "      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > " << drmax_bins[drmax_bin] << " && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < " << drmax_bins[drmax_bin+1] << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > " << drmax_bins[drmax_bin] << " && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < " << drmax_bins[drmax_bin+1] << ") \n";
    }
    comp_output_file << "        weight = weight*" << zll_variations[1][drmax_bin] <<";\n";
  }
  first_bin = true;
  for (unsigned int met_bin = 0; met_bin < met_bins.size()-1; met_bin++) {
    //set the upper bound to the edge of the bin or 9999 if top bin
    double upper_bound = met_bins[met_bin+1];
    if (met_bin == met_bins.size()-2) {
      upper_bound = 9999;
    }
    if (first_bin) {
      comp_output_file << "      if (b.met() > " << met_bins[met_bin] << " && b.met() < " << upper_bound << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (b.met() > " << met_bins[met_bin] << " && b.met() < " << upper_bound << ") \n";
    }
    comp_output_file << "        weight = weight*" << zll_variations[2][met_bin] <<";\n";
  }
  comp_output_file << "      break;\n";
  comp_output_file << "    case 7:\n";
  first_bin = true;
  for (unsigned int nb_bin = 0; nb_bin < nb_bins.size()-1; nb_bin++) {
    if (first_bin) {
      comp_output_file << "      if (Higfuncs::jetid_nb.GetScalar(b) > " << nb_bins[nb_bin] << " && Higfuncs::jetid_nb.GetScalar(b) < " << nb_bins[nb_bin+1] << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (Higfuncs::jetid_nb.GetScalar(b) > " << nb_bins[nb_bin] << " && Higfuncs::jetid_nb.GetScalar(b) < " << nb_bins[nb_bin+1] << ") \n";
    }
    comp_output_file << "        weight = weight*" << qcd_variations[0][nb_bin] <<";\n";
  }
  first_bin = true;
  for (unsigned int drmax_bin = 0; drmax_bin < drmax_bins.size()-1; drmax_bin++) {
    if (first_bin) {
      comp_output_file << "      if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > " << drmax_bins[drmax_bin] << " && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < " << drmax_bins[drmax_bin+1] << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (Higfuncs::jetid_hig_cand_drmax.GetScalar(b) > " << drmax_bins[drmax_bin] << " && Higfuncs::jetid_hig_cand_drmax.GetScalar(b) < " << drmax_bins[drmax_bin+1] << ") \n";
    }
    comp_output_file << "        weight = weight*" << qcd_variations[1][drmax_bin] <<";\n";
  }
  first_bin = true;
  for (unsigned int met_bin = 0; met_bin < met_bins.size()-1; met_bin++) {
    //set the upper bound to the edge of the bin or 9999 if top bin
    double upper_bound = met_bins[met_bin+1];
    if (met_bin == met_bins.size()-2) {
      upper_bound = 9999;
    }
    if (first_bin) {
      comp_output_file << "      if (b.met() > " << met_bins[met_bin] << " && b.met() < " << upper_bound << ") \n";
      first_bin = false;
    }
    else {
      comp_output_file << "      else if (b.met() > " << met_bins[met_bin] << " && b.met() < " << upper_bound << ") \n";
    }
    comp_output_file << "        weight = weight*" << qcd_variations[2][met_bin] <<";\n";
  }
  comp_output_file << "    default:\n";
  comp_output_file << "      break;\n";
  comp_output_file << "  }\n";
  comp_output_file << "  return weight;\n";
  comp_output_file << "});\n";
  comp_output_file.close();
  std::cout << "Wrote src/higgsino/" << options.out_filename << std::endl;

  time(&endtime); 
  std::cout << std::endl << "Took " << difftime(endtime, begtime) << " seconds\n" << std::endl;
}

std::vector<std::vector<double>> get_variations_from_ratio(PlotMaker &pm, unsigned int pm_index) {
  if (pm.Figures().size() < pm_index+2) {
    std::cout << "ERROR: insufficient plots in plot maker" << std::endl;
    return {};
  }
  double bot_min = 0, bot_max = 0;
  //nb variations
  Hist1D *nb_variations_hist = static_cast<Hist1D*>(pm.Figures()[pm_index].get());
  TH1D nb_variations_th1 = nb_variations_hist->GetBottomPlots(bot_min, bot_max).at(0);
  std::vector<double> nb_variations;
  for (unsigned int nb_bin = 1; nb_bin < nb_bins.size(); nb_bin++) {
    nb_variations.push_back(nb_variations_th1.GetBinContent(static_cast<int>(nb_bin)));
  }
  //drmax variations
  Hist1D *drmax_variations_hist = static_cast<Hist1D*>(pm.Figures()[pm_index+1].get());
  TH1D drmax_variations_th1 = drmax_variations_hist->GetBottomPlots(bot_min, bot_max).at(0);
  std::vector<double> drmax_variations;
  for (unsigned int drmax_bin = 1; drmax_bin < drmax_bins.size(); drmax_bin++) {
    drmax_variations.push_back(drmax_variations_th1.GetBinContent(static_cast<int>(drmax_bin)));
  }
  //met variations
  Hist1D *met_variations_hist = static_cast<Hist1D*>(pm.Figures()[pm_index+2].get());
  TH1D met_variations_th1 = met_variations_hist->GetBottomPlots(bot_min, bot_max).at(0);
  std::vector<double> met_variations;
  for (unsigned int met_bin = 1; met_bin < met_bins.size(); met_bin++) {
    met_variations.push_back(met_variations_th1.GetBinContent(static_cast<int>(met_bin)));
  }
  return {nb_variations, drmax_variations, met_variations};
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

