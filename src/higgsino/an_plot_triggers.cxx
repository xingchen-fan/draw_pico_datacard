//This script generates the plots seen in the trigger [systematics] section of 
//the AN as well as the tables of trigger systematics and C++ files that can be
//used to apply trigger efficiencies
//
//Arguments
// --single_thread (-s) to run single thread for debugging
// --unblind (-u) to unblind (not done for AN plots)
// --year (-y) yearname to select data year (2016, 2017, 2018, run2)
// --tag (-t) to add a tag to produced plots
// --string_options (-o) to specify what to do
//
//possible string options (default indicates it is an AN plot): 
// plot (default) - generate AN plots except extrapolation plots
// systematic - generate systematics table
// efficiency - generate trigger efficiency files
// cr - also generate plots and efficiencies for control regions (much slower)

#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/efficiency_plot.hpp"
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"
#include "higgsino/hig_utilities.hpp"
#include "higgsino/script_utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  //constants
  //bins for systematics
  const std::vector<double> sys_met_bins = {150,160,180,200,225,250,300,350};
  const std::vector<double> sys_ht_bins = {0, 300, 400, 600, 800, 9999};
  //bins for search region 
  const std::vector<double> search_ht_bins = {0,200,300,400,600,950,9999};
  const std::vector<std::vector<double>> search_met_bins = {
      {150, 155, 160, 165, 170, 180, 190, 200, 350}, //HT 0-200
      {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 350}, //HT 200-300
      {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 350}, //HT 300-400
      {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 350}, //HT 400-600
      {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 350}, //HT 600-950
      {150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 275, 300, 350}}; //HT 950+
  //bins for fake met control region
  const std::vector<double> qcd_ht_bins = {0,350,450,550,650,800,1000,9999};
  const std::vector<std::vector<double>> qcd_met_bins = {
      {150,160,170,180,190,200,225,250,600}, //HT 0-350
      {150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,300,600}, //HT 350-450
      {150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,300,400,600}, //HT 450-550
      {150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,400,600}, //HT 550-650
      {150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350,400,450,500,600}, //HT 650-800
      {150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350,400,450,500,600}, //HT 800-1000
      {150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350,400,450,500,600}}; //HT 1000+
  //bins for 1l control region
  const std::vector<double> singlelep_met_bins{0,50,75,100,125,150,175,215,300};
  const std::vector<double> singlelep_ht_bins{0,400,600,800};
  const std::vector<double> el_pt_bins{20,25,30,40,110,120,150};
  const std::vector<double> mu_pt_bins{20,25,30,50,100};
  //bins for 2l control region
  const std::vector<double> twoel_pt_bins{40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,120};
  const std::vector<double> twomu_pt_bins{40,45,50,60};
}

struct EfficiencyVariable {
  std::vector<std::vector<double>> variable_bins;
  std::vector<double> bins;
  std::string name;
  std::string expression;
};

//systematic helper function declarations
std::vector<double> trig_postprocess_systematics(PlotMaker* pm, script_utilities::ArgStruct options);
std::vector<std::vector<double>> generate_syst_extrapolation_table(PlotMaker* pm, std::shared_ptr<Process> proc, script_utilities::ArgStruct options);
std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name);
std::vector<double> get_systematic_table_values(TGraphAsymmErrors *gr);
std::vector<double> get_systematic_table_variations(std::vector<std::vector<double>> sys_table);
std::string make_systematic_table_row(std::vector<double> values, bool total=false);
std::vector<double> add_variations_quadrature(std::vector<std::vector<double>> values);
std::vector<double> max_variations(std::vector<std::vector<double>> values);

//efficiency helper function declarations
void trig_postprocess_efficiencies(PlotMaker* pm, 
                                   script_utilities::ArgStruct options,
                                   std::vector<double> systematics, 
                                   std::vector<std::vector<double>> extrapolation_systematics);
void generate_1d_efficiency(PlotMaker* pm, std::string func_name, std::string plot_name,
                            EfficiencyVariable x_var, ofstream &out_file);
void generate_2d_efficiency_variable_binning(PlotMaker* pm, 
                                             std::string func_name,
                                             std::vector<std::string> plot_name,
                                             EfficiencyVariable x_var, 
                                             EfficiencyVariable y_var,
                                             ofstream &out_file, 
                                             std::vector<double> systematics, 
                                             std::vector<std::vector<double>> extrapolation_systematics);
void generate_3d_efficiency(PlotMaker* pm, std::string func_name,
                            std::vector<std::string> plot_name,
                            EfficiencyVariable x_var, EfficiencyVariable y_var, 
                            EfficiencyVariable z_var, ofstream &out_file);

int main(int argc, char *argv[]){
  //------------------------------------------------------------------------------------
  //                                    initialization
  //------------------------------------------------------------------------------------
  gErrorIgnoreLevel = 6000;
  time_t begtime, endtime;
  time(&begtime);

  script_utilities::ArgStruct options = script_utilities::get_options(
      argc, argv, "plot,cr");
  
  std::string out_filename = "triggereff";

  // Define 1D+2D plot types of interest; include overflows for efficiencies
  Palette colors("txt/colors.txt", "default");
  PlotOpt lin_lumi("txt/plot_styles.txt", "Std1D");
  vector<PlotOpt> all_plot_types = {lin_lumi.Title(TitleType::info).Overflow(OverflowType::overflow)};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).Overflow(OverflowType::overflow)};
  //samples
  std::set<int> years;
  HigUtilities::parseYears(options.year_string, years);
  string total_luminosity_string = HigUtilities::getLuminosityString(options.year_string);
  if (HigUtilities::is_in_string_options(options.string_options,"systematic")
      || HigUtilities::is_in_string_options(options.string_options,"efficiency")) {
    if (years.size() > 1) {
      ERROR("Can only have 1 year when using the systematic or efficiency options.");
      return 1;
    }
  }

  string base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_klamath/";
  string data_skim_met150_folder = "data/skim_met150/";
  string data_skim_1l2j_folder = "data/skim_1l2j/";
  string mc_skim_met150_folder = "mc/skim_met150/";

  vector<shared_ptr<Process> > procs_met150 = {
      Process::MakeShared<Baby_pico>("Data", Process::Type::data, 
      kBlack, attach_folder(base_folder, years, data_skim_met150_folder, 
      {"*.root"}), "1")};
  vector<shared_ptr<Process> > procs_met150_mc = {
      Process::MakeShared<Baby_pico>("t#bar{t} 1l", Process::Type::background, 
      kBlack, attach_folder(base_folder, years, mc_skim_met150_folder, 
      {"*TTJets_SingleLeptFromT*.root"}), "stitch")};
  vector<shared_ptr<Process> > procs_1l2j = {
      Process::MakeShared<Baby_pico>("Data", Process::Type::data, 
      kBlack, attach_folder(base_folder, years, data_skim_1l2j_folder, 
      {"*.root"}), "1")};

  // Set MC 
  std::map<std::string, std::set<std::string>> mctags; 
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
  mctags["all"] = std::set<std::string>({"*TTJets_SingleLept*",
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

  std::vector<std::shared_ptr<Process>> procs_data_met150_allyears = {};
  std::vector<std::shared_ptr<Process>> procs_mc_met150_allyears = {};
  std::vector<std::shared_ptr<Process>> procs_data_1l2j_allyears = {};
  std::vector<std::vector<std::shared_ptr<Process>>> procs_data_ambb = {};
  std::vector<std::string> year_strings;
  std::vector<std::string> luminosity_year_string;

  unsigned int year_idx = 0;
  for (int iyear : years) {
    std::string iyear_string = "2018";
    std::string iyear_lumi = "60";
    std::set<int> iyear_set = {iyear};
    short iyear_color = kRed;
    if (iyear==2016) { 
      iyear_string = "2016";
      iyear_lumi = "35.9";
      iyear_color = kBlue;
    }
    else if (iyear==2017) {
      iyear_string = "2017";
      iyear_lumi = "41.5";
      iyear_color = kGreen;
    }
    procs_data_met150_allyears.push_back(
        Process::MakeShared<Baby_pico>("Data "+iyear_string, Process::Type::data, iyear_color,
        attach_folder(base_folder, iyear_set, data_skim_met150_folder, {"*.root"}), "stitch"));
    procs_mc_met150_allyears.push_back(
        Process::MakeShared<Baby_pico>("Data "+iyear_string, Process::Type::data, iyear_color,
        attach_folder(base_folder, iyear_set, mc_skim_met150_folder, mctags["all"]), "stitch"));
    procs_data_1l2j_allyears.push_back(
        Process::MakeShared<Baby_pico>("Data "+iyear_string, Process::Type::data, iyear_color,
        attach_folder(base_folder, iyear_set, data_skim_1l2j_folder, {"*.root"}), "stitch"));
    procs_data_ambb.push_back(std::vector<std::shared_ptr<Process>>());
    procs_data_ambb[year_idx].push_back(
        Process::MakeShared<Baby_pico>("Data "+iyear_string+" #LT m_{bb} #GT < 100 GeV", 
        Process::Type::data, kRed,
        attach_folder(base_folder, iyear_set, data_skim_met150_folder, {"*.root"}), 
        "stitch&&njet>=4&&hig_cand_am[0]<100"));
    procs_data_ambb[year_idx].push_back(
        Process::MakeShared<Baby_pico>("Data "+iyear_string+" 100 #leq #LT m_{bb} #GT < 150 GeV", 
        Process::Type::data, kGreen,
        attach_folder(base_folder, iyear_set, data_skim_met150_folder, {"*.root"}), 
        "stitch&&njet>=4&&hig_cand_am[0]>=100&&hig_cand_am[0]<150"));
    procs_data_ambb[year_idx].push_back(
        Process::MakeShared<Baby_pico>("Data "+iyear_string+"#LT m_{bb} #GT #geq 150 GeV", 
        Process::Type::data, kBlue,
        attach_folder(base_folder, iyear_set, data_skim_met150_folder, {"*.root"}), 
        "stitch&&njet>=4&&hig_cand_am[0]>=150"));
    year_strings.push_back(iyear_string);
    luminosity_year_string.push_back(iyear_lumi);
    year_idx += 1;
  }

  //------------------------------------------------------------------------------------
  //                                     named funcs
  //------------------------------------------------------------------------------------

  const NamedFunc high_pt_jet("high_pt_jet", [](const Baby &b) -> NamedFunc::ScalarType{
    bool r_high_pt_jet = false;
    if (b.njet()>0) {
      if (b.jet_pt()->at(0) > 500) {
        r_high_pt_jet = true;
      }
    }
    return r_high_pt_jet;
  });

  const NamedFunc met("met", [](const Baby &b) -> NamedFunc::ScalarType{
      return b.met();
  });

  const NamedFunc ht("ht", [](const Baby &b) -> NamedFunc::ScalarType{
      return b.ht();
  });

  const NamedFunc nb_higgsino("nb_higgsino", [](const Baby &b) -> NamedFunc::ScalarType{
    int r_nb_higgsino;
    if (b.nbm() == 0) r_nb_higgsino = 0;
    else if (b.nbm() == 1) r_nb_higgsino = 1;
    else if (b.nbt() == 2 && b.nbm() == 2) r_nb_higgsino = 2;
    else if (b.nbt() >= 2 && b.nbm() == 3 && b.nbl() == 3) r_nb_higgsino = 3;
    else if (b.nbt() >= 2 && b.nbm() >= 3 && b.nbl() >= 4) r_nb_higgsino = 4;
    else r_nb_higgsino = -1;
    return r_nb_higgsino;
  });

  //------------------------------------------------------------------------------------
  //                                          plots
  //------------------------------------------------------------------------------------

  PlotMaker pm;

  //plots of efficiency with respect to various analysis variables
  if (HigUtilities::is_in_string_options(options.string_options,"plot")) {
    //-------------6.2 table 1 plots (0l Real MET: MET, HT(low MET), HT(high MET))-------------
    pm.Push<EfficiencyPlot>(Axis(100, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && 
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_met").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&150<met&&met<=200" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_ht_lowmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met&&met<=300" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_ht_highmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //-------------6.3 table 1 plots (0l Fake MET: MET, HT(low MET), HT(high MET))-------------
    pm.Push<EfficiencyPlot>(Axis(100, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger,
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lfakemet_met").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
        Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&150<met&&met<=200" && Higfuncs::jetht_trigger,
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lfakmet_ht_lowmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
        Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&200<met&&met<=300" && Higfuncs::jetht_trigger,
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lfakmet_ht_highmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //-------------9.3 table 1 plots (0l Real MET: Nj Nb DRmax, Ambb(not in AN))-------------
    pm.Push<EfficiencyPlot>(Axis(5, -0.5, 4.5, nb_higgsino, "N_{b}", {}),
        Higfuncs::pass_filters && "njet>=4&&njet<=5&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_nb").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(5, 2.5, 7.5, "njet", "N_{j}", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_njet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(20, 0.0, 4.0, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
        Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_drmax").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    pm.Push<EfficiencyPlot>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {}),
        Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig__eff_0lrealmet_ambb");
    //-------------9.3 table 2 plots (0l Real MET: <mbb> dependence)-------------
    for (unsigned int year_index = 0; year_index < procs_data_ambb.size(); year_index++) {
      pm.Push<EfficiencyPlot>(Axis(30, 150., 300., "met", "Offline p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300" &&
          (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
          Higfuncs::met_trigger,
          procs_data_ambb[year_index]).Tag("FixName:trig__eff_0lrealmet_met_ht250to300_ambbcuts_"+year_strings[year_index]).YTitle("p_{T}^{miss} Triggers").LuminosityTag(luminosity_year_string[year_index]);
      pm.Push<EfficiencyPlot>(Axis(30, 150., 300., "met", "Offline p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400" &&
          (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
          Higfuncs::met_trigger,
          procs_data_ambb[year_index]).Tag("FixName:trig__eff_0lrealmet_met_ht300to400_ambbcuts_"+year_strings[year_index]).YTitle("p_{T}^{miss} Triggers").LuminosityTag(luminosity_year_string[year_index]);
      pm.Push<EfficiencyPlot>(Axis(30, 150., 300., "met", "Offline p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600" &&
          (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
          Higfuncs::met_trigger,
          procs_data_ambb[year_index]).Tag("FixName:trig__eff_0lrealmet_met_ht400to600_ambbcuts_"+year_strings[year_index]).YTitle("p_{T}^{miss} Triggers").LuminosityTag(luminosity_year_string[year_index]);
    }
    if (HigUtilities::is_in_string_options(options.string_options,"cr")) {
      //-------------6.4 table 1 plots (1l CR: MET lep_pt HT(lowmet) HT(highmet))-------------
      //electrons
      pm.Push<EfficiencyPlot>(Axis(55, 0., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==1" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1el_met").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Electron p_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1el_leppt").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1el_ht_lowmet").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==1&&200<=met&&met<300" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1el_ht_highmet").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
      //muons
      pm.Push<EfficiencyPlot>(Axis(55, 0., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==1" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1mu_met").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Muon p_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1mu_leppt").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1mu_ht_lowmet").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==1&&200<=met&&met<300" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__1mu_ht_highmet").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
      //-------------6.5 table 1 plots (2l CR: lep_pt HT)-------------
      //electrons
      pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::el_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__2el_leppt").YTitle("Electron Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::el_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__2el_ht").YTitle("Electron Triggers").LuminosityTag(total_luminosity_string);
      //muons
      pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::mu_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__2mu_leppt").YTitle("Muon Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::mu_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig__eff__2mu_ht").YTitle("Muon Triggers").LuminosityTag(total_luminosity_string);
    }
    ////-------------Data vs MC comparison, not in AN-------------
    //pm.Push<EfficiencyPlot>(Axis(100, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht<250" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    "HLT_PFMET120_PFMHT120_IDTight",
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_onlymet120").YTitle("HLT_MET120_MHT120_IDTight").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(100, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht<250" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    "HLT_PFMET120_PFMHT120_IDTight",
    //    procs_mc_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_mc_onlymet120").YTitle("HLT_MET120_MHT120_IDTight").LuminosityTag(total_luminosity_string);
    //-------------HT binned plots to determine trigger binning for 0lrealmet-------------
    //pm.Push<EfficiencyPlot>(Axis(50, 150., 400., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht<200" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_ht0to200").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(50, 150., 400., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht>=200&&ht<300" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_ht200to300").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(50, 150., 400., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht>=300&&ht<400" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_ht300to400").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(50, 150., 400., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht>=400&&ht<600" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_ht400to600").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(50, 150., 400., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht>=600&&ht<950" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_ht600to950").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(50, 150., 400., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht>=950" && Higfuncs::el_trigger && (Higfuncs::lead_signal_lepton_pt < 35), 
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met_ht950toInf").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    ////-------------HT binned plots to determine trigger binning for 0lfakemet-------------
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=0&&ht<350",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht0to350").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=350&&ht<450",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht350to450").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=450&&ht<550",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht450to550").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=550&&ht<650",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht550to650").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=650&&ht<800",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht650to800").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=800&&ht<1000",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht800to1000").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //pm.Push<EfficiencyPlot>(Axis(80, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
    //    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger && "ht>=1000",
    //    Higfuncs::met_trigger,
    //    procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met_ht1000toInf").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  }
  //0l systematics plots
  if (HigUtilities::is_in_string_options(options.string_options,"systematic")) {
    //nj3 (nominal)
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_nj3_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //nj4
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet==4&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_nj4_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //nj5
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet==5&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_nj5_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //nb0
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&nbt==0" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_nb0_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //nb1
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&nbt==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_nb1_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //nb2
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&nbt==2&&nbm==2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_nb2_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //drmax<2.2
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]<2.2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_lowdrmax_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //drmax>2.2
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]>=2.2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_highdrmax_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //300<=ht<400 (nominal)
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_ht300to400_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //300<=ht<400, am<140
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]<140" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_ht300to400_lowam_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //300<=ht<400, am>140
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]>=140" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_ht300to400_higham_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //jet1 pt>500 (nominal)
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && high_pt_jet, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_highptjet_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //jet1 pt>500 ref trigger jet500
    pm.Push<EfficiencyPlot>(Axis(sys_met_bins, "met", "p_{T}^{miss} [GeV]", {}),
        Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&HLT_PFJet500" && (Higfuncs::lead_signal_lepton_pt<35) && high_pt_jet, 
        Higfuncs::met_trigger,
        procs_met150).Tag("FixName:trig__sys_highptjet_jetref_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    //extrapolation systematic
    for (unsigned int ht_bin_idx = 0; ht_bin_idx < (sys_ht_bins.size()-1); ht_bin_idx++) {
      for (unsigned int met_bin_idx = 0; met_bin_idx < (sys_met_bins.size()-1); met_bin_idx++) {
        pm.Push<EfficiencyPlot>(Axis(16, 0., 80., Higfuncs::lead_signal_lepton_pt, "Electron p_{T} [GeV]", {}),
            Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && Higfuncs::el_trigger && 
            (ht >= sys_ht_bins[ht_bin_idx]) && (ht < sys_ht_bins[ht_bin_idx+1]) && 
            (met >= sys_met_bins[met_bin_idx]) && (met < sys_met_bins[met_bin_idx+1]),
            Higfuncs::met_trigger,
            procs_met150).Tag("FixName:trig__sys_0l_htbin"+std::to_string(ht_bin_idx)+"_metbin"+std::to_string(met_bin_idx)+"_"+options.year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
        //TODO: fix to allow generating efficiencies for multiple years at once?
      }
    }
  }

  // trigger efficiency plots and application file
  if (HigUtilities::is_in_string_options(options.string_options,"efficiency")) {
    //0l Real MET 
    for (unsigned int ht_bin_idx = 0; ht_bin_idx < (search_ht_bins.size()-1); ht_bin_idx++) {
      std::vector<double> met_bins = search_met_bins[ht_bin_idx];
      pm.Push<EfficiencyPlot>(Axis(met_bins, "met", "p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && 
          (Higfuncs::lead_signal_lepton_pt<35)  && Higfuncs::el_trigger && 
          (ht >= search_ht_bins[ht_bin_idx]) && (ht < search_ht_bins[ht_bin_idx+1]),
          Higfuncs::met_trigger, procs_met150)
          .Tag("FixName:trig__eff_search_htbin"+std::to_string(ht_bin_idx)+"_"+options.year_string)
          .YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    }
    //0l Fake MET
    for (unsigned int ht_bin_idx = 0; ht_bin_idx < (qcd_ht_bins.size()-1); ht_bin_idx++) {
      std::vector<double> met_bins = qcd_met_bins[ht_bin_idx];
      pm.Push<EfficiencyPlot>(Axis(met_bins, "met", "p_{T}^{miss} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && 
          (Higfuncs::lead_signal_lepton_pt<35)  && Higfuncs::el_trigger && 
          (ht >= qcd_ht_bins[ht_bin_idx]) && (ht < qcd_ht_bins[ht_bin_idx+1]),
          Higfuncs::met_trigger, procs_met150)
          .Tag("FixName:trig__eff_qcd_htbin"+std::to_string(ht_bin_idx)+"_"+options.year_string)
          .YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
    }
    if (HigUtilities::is_in_string_options(options.string_options,"cr")) {
      //Single electron
      for (unsigned int ht_bin_idx = 0; ht_bin_idx < (singlelep_ht_bins.size()-1); ht_bin_idx++) {
        for (unsigned int met_bin_idx = 0; met_bin_idx < (singlelep_met_bins.size()-1); met_bin_idx++) {
          pm.Push<EfficiencyPlot>(Axis(el_pt_bins, Higfuncs::lead_signal_lepton_pt, "electron p_{T} [GeV]", {}),
              Higfuncs::pass_filters && "njet>=2&&nel==1" && Higfuncs::jetht_trigger &&
              (ht >= singlelep_ht_bins[ht_bin_idx]) && (ht < singlelep_ht_bins[ht_bin_idx+1]) &&
              (met >= singlelep_met_bins[met_bin_idx]) && (met < singlelep_met_bins[met_bin_idx+1]),
              Higfuncs::met_trigger || Higfuncs::el_trigger, procs_1l2j)
              .Tag("FixName:trig__eff_1e_htbin"+std::to_string(ht_bin_idx)+
              "_metbin"+std::to_string(met_bin_idx)+"_"+options.year_string)
              .YTitle("p_{T}^{miss} OR Single-e Triggers").LuminosityTag(total_luminosity_string);
        }
      }
      //Single muon
      for (unsigned int ht_bin_idx = 0; ht_bin_idx < (singlelep_ht_bins.size()-1); ht_bin_idx++) {
        for (unsigned int met_bin_idx = 0; met_bin_idx < (singlelep_met_bins.size()-1); met_bin_idx++) {
          pm.Push<EfficiencyPlot>(Axis(mu_pt_bins, Higfuncs::lead_signal_lepton_pt, "muon p_{T} [GeV]", {}),
              Higfuncs::pass_filters && "njet>=2&&nmu==1" && Higfuncs::jetht_trigger &&
              (ht >= singlelep_ht_bins[ht_bin_idx]) && (ht < singlelep_ht_bins[ht_bin_idx+1]) &&
              (met >= singlelep_met_bins[met_bin_idx]) && (met < singlelep_met_bins[met_bin_idx+1]),
              Higfuncs::met_trigger || Higfuncs::mu_trigger, procs_1l2j)
              .Tag("FixName:trig__eff_1mu_htbin"+std::to_string(ht_bin_idx)+
              "_metbin"+std::to_string(met_bin_idx)+"_"+options.year_string)
              .YTitle("p_{T}^{miss} OR Single-#mu Triggers").LuminosityTag(total_luminosity_string);
        }
      }
      //Two leptons
      pm.Push<EfficiencyPlot>(Axis(twoel_pt_bins, Higfuncs::lead_signal_lepton_pt, "leading electron p_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && 
          (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::el_trigger, procs_1l2j)
          .Tag("FixName:trig__eff_2e_"+options.year_string)
          .YTitle("Single-e Triggers").LuminosityTag(total_luminosity_string);
      pm.Push<EfficiencyPlot>(Axis(twomu_pt_bins, Higfuncs::lead_signal_lepton_pt, "leading muon p_{T} [GeV]", {}),
          Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && 
          (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::mu_trigger, procs_1l2j)
          .Tag("FixName:trig__eff_2mu_"+options.year_string)
          .YTitle("Single-#mu Triggers").LuminosityTag(total_luminosity_string);
    }
  }

  pm.multithreaded_ = !options.single_thread;
  pm.min_print_ = true;
  pm.print_2d_figures_ = false;
  pm.MakePlots(1.0);

  //-------------------------------------------------------------------
  // post processing
  //-------------------------------------------------------------------

  //also save plots to root file in case something fails
  //TFile* out_file = TFile::Open(("ntuples/"+out_filename+".root").c_str(),"RECREATE");
  std::vector<double> systematics;
  std::vector<std::vector<double>> extrapolation_systematics;
  if (HigUtilities::is_in_string_options(options.string_options,"systematic")) {
    systematics = trig_postprocess_systematics(&pm, options);
    extrapolation_systematics = generate_syst_extrapolation_table(&pm, procs_met150[0], options);
  }
  if (HigUtilities::is_in_string_options(options.string_options,"efficiency")) {
    std::cout << "DEBUG: starting efficiency" << std::endl;
    trig_postprocess_efficiencies(&pm, options, systematics, extrapolation_systematics);
    std::cout << "DEBUG: finishing efficiency" << std::endl;
  }
  //out_file->Close();

  time(&endtime);
  cout<<endl<<"Processing took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

//---------------------------------------------------------------
// systematics functions
//---------------------------------------------------------------
std::vector<double> trig_postprocess_systematics(PlotMaker* pm, script_utilities::ArgStruct options) {
  //generate ratio plots
  EfficiencyPlot * eff_nj3 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_nj3_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_nj3 = static_cast<TGraphAsymmErrors*>((eff_nj3->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_nj4 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_nj4_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_nj4 = static_cast<TGraphAsymmErrors*>((eff_nj4->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_nj5 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_nj5_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_nj5 = static_cast<TGraphAsymmErrors*>((eff_nj5->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_nb0 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_nb0_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_nb0 = static_cast<TGraphAsymmErrors*>((eff_nb0->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_nb1 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_nb1_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_nb1 = static_cast<TGraphAsymmErrors*>((eff_nb1->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_nb2 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_nb2_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_nb2 = static_cast<TGraphAsymmErrors*>((eff_nb2->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_lowdrmax = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_lowdrmax_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_lowdrmax = static_cast<TGraphAsymmErrors*>((eff_lowdrmax->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_highdrmax = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_highdrmax_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_highdrmax = static_cast<TGraphAsymmErrors*>((eff_highdrmax->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_ht300to400 = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_ht300to400_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_ht300to400 = static_cast<TGraphAsymmErrors*>((eff_ht300to400->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_ht300to400_lowam = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_ht300to400_lowam_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_ht300to400_lowam = static_cast<TGraphAsymmErrors*>((eff_ht300to400_lowam->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_ht300to400_higham = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_ht300to400_higham_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_ht300to400_higham = static_cast<TGraphAsymmErrors*>((eff_ht300to400_higham->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_highptjet = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_highptjet_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_highptjet = static_cast<TGraphAsymmErrors*>((eff_highptjet->data_ratio_plots_[0]).get()->Clone());
  EfficiencyPlot * eff_highptjet_jetref = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_highptjet_jetref_"+options.year_string).get());
  TGraphAsymmErrors* plot_ratio_highptjet_jetref = static_cast<TGraphAsymmErrors*>((eff_highptjet_jetref->data_ratio_plots_[0]).get()->Clone());
  
  //bool do_datamc = false;
  ofstream table_file;
  std::vector<std::vector<double>> variations_table;
  std::vector<double> variations_row;
  std::vector<std::vector<double>> ambb_table;
  std::vector<double> ambb_row;
  std::vector<std::vector<double>> reference_table;
  std::vector<double> reference_row;
  std::vector<std::vector<double>> nel_table;
  std::vector<double> nel_row;
  std::vector<std::vector<double>> datamc_table;
  std::vector<double> datamc_row;
  std::vector<std::vector<double>> greater_variations_rows;
  std::vector<std::vector<double>> systematics_rows;
  table_file.open("tables/trig__sys_"+options.year_string+".tex");
  table_file << "\\begin{tabular}{lccccccc}\\hline \\hline \n";
  table_file << "p$_\\text{T}^\\text{miss}$ range& [150,160]& [160,180]& [180,200]& [200,225]& [225,250]& [250,300]& [300+]\\\\ \\hline \\hline \n";
  table_file << "\\multicolumn{8}{c}{Kinematic variations}\\\\ \\hline \n";
  table_file << make_systematic_table_row(plot_ratio_nj3,"$N_\\text{jets}\\geq 3$ (Nominal)") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_nj3));
  table_file << make_systematic_table_row(plot_ratio_nj4,"$N_\\text{jets}=4$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_nj4));
  table_file << make_systematic_table_row(plot_ratio_nj5,"$N_\\text{jets}=5$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_nj5));
  table_file << make_systematic_table_row(plot_ratio_nb0,"$N_\\text{b}=0$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_nb0));
  table_file << make_systematic_table_row(plot_ratio_nb1,"$N_\\text{b}=1$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_nb1));
  table_file << make_systematic_table_row(plot_ratio_nb2,"$N_\\text{b}\\geq2$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_nb2));
  table_file << make_systematic_table_row(plot_ratio_lowdrmax,"$\\Delta R_\\text{max}< 2.2$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_lowdrmax));
  table_file << make_systematic_table_row(plot_ratio_highdrmax,"$\\Delta R_\\text{max}\\geq 2.2$") << " \\hline \n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_highdrmax));
  variations_row = get_systematic_table_variations(variations_table);
  greater_variations_rows.push_back(variations_row);
  table_file << make_systematic_table_row(variations_row) << "\n";
  table_file << "\\multicolumn{8}{c}{$\\langle m_{bb}\\rangle$ variations (for $300 \\leq H_{T} < 400$~GeV)}\\\\ \\hline \n";
  table_file << make_systematic_table_row(plot_ratio_ht300to400,"$\\langle m_{bb}\\rangle$ Inclusive (Nominal)") << "\n";
  ambb_table.push_back(get_systematic_table_values(plot_ratio_ht300to400));
  table_file << make_systematic_table_row(plot_ratio_ht300to400_lowam,"$\\langle m_{bb}\\rangle \\leq$ 140 GeV") << "\n";
  ambb_table.push_back(get_systematic_table_values(plot_ratio_ht300to400_lowam));
  table_file << make_systematic_table_row(plot_ratio_ht300to400_higham,"$\\langle m_{bb}\\rangle >$ 140 GeV") << " \\hline \n";
  ambb_table.push_back(get_systematic_table_values(plot_ratio_ht300to400_higham));
  ambb_row = get_systematic_table_variations(ambb_table);
  greater_variations_rows.push_back(ambb_row);
  table_file << make_systematic_table_row(ambb_row) << "\n";
  systematics_rows.push_back(max_variations(greater_variations_rows));
  table_file << make_systematic_table_row(max_variations(greater_variations_rows)) << "\n";
  table_file << "\\multicolumn{8}{c}{Reference trigger (for Max $p_\\text{Tjet}>500$~GeV)}\\\\ \\hline \n";
  table_file << make_systematic_table_row(plot_ratio_highptjet,"SingleElectron (Nominal)") << "\n";
  reference_table.push_back(get_systematic_table_values(plot_ratio_highptjet));
  table_file << make_systematic_table_row(plot_ratio_highptjet_jetref,"Jet500") << " \\hline \n";
  reference_table.push_back(get_systematic_table_values(plot_ratio_highptjet_jetref));
  reference_row = get_systematic_table_variations(reference_table);
  systematics_rows.push_back(reference_row);
  table_file << make_systematic_table_row(reference_row) << "\n";
  //if (do_datamc) {
  //  table_file << "\\multicolumn{8}{c}{Data compared to MC (MET120\\_HT60 trigger only)}\\\\ \\hline \n";
  //  table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio")),"Data (nominal)") << "\n";
  //  datamc_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio"))));
  //  table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio")),"MC") << "\\hline \n";
  //  datamc_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio"))));
  //  datamc_row = get_systematic_table_variations(datamc_table);
  //  systematics_rows.push_back(datamc_row);
  //  table_file << make_systematic_table_row(datamc_row) << "\n";
  //}
  std::vector<double> ret = add_variations_quadrature(systematics_rows);
  table_file << make_systematic_table_row(ret, true) << "\n";
  table_file << "\\end{tabular} \n";
  table_file.close();
  std::cout << "Wrote tables/higgsino/trig__sys_"+options.year_string+".tex" << std::endl;
  delete plot_ratio_nj3;
  delete plot_ratio_nj4;
  delete plot_ratio_nj5;
  delete plot_ratio_nb0;
  delete plot_ratio_nb1;
  delete plot_ratio_nb2;
  delete plot_ratio_lowdrmax;
  delete plot_ratio_highdrmax;
  delete plot_ratio_ht300to400;
  delete plot_ratio_ht300to400_lowam;
  delete plot_ratio_ht300to400_higham;
  delete plot_ratio_highptjet;
  delete plot_ratio_highptjet_jetref;
  return ret;
}

std::vector<std::vector<double>> generate_syst_extrapolation_table(PlotMaker* pm, std::shared_ptr<Process> proc, script_utilities::ArgStruct options) {
  //generate extrapolation systematics
  std::cout << "DEBUG: starting extrapolation systematics" << std::endl;
  TCanvas * c = new TCanvas("c");
  TF1 linear = TF1("linear","[0]*x+[1]",0,50,"");
  std::vector<std::vector<double>> nominal_values;
  std::vector<std::vector<double>> extrapolated_values;
  std::vector<std::vector<double>> extrapolation_syst;
  for (unsigned int ht_bin_idx = 0; ht_bin_idx < (sys_ht_bins.size()-1); ht_bin_idx++) {
    nominal_values.push_back({});
    extrapolated_values.push_back({});
    extrapolation_syst.push_back({});
    for (unsigned int met_bin_idx = 0; met_bin_idx < (sys_met_bins.size()-1); met_bin_idx++) {
      std::cout << "DEBUG: extrapolation MET bin " << met_bin_idx << std::endl;
      EfficiencyPlot * eff_plot = static_cast<EfficiencyPlot*>(pm->GetFigure("trig__sys_0l_htbin"+std::to_string(ht_bin_idx)+"_metbin"+std::to_string(met_bin_idx)+"_"+options.year_string).get());
      EfficiencyPlot::SingleEfficiencyPlot* singleff_plot = static_cast<EfficiencyPlot::SingleEfficiencyPlot*>(eff_plot->GetComponent(proc.get()));
      TH1D den_hist = singleff_plot->raw_denominator_hist_;
      TH1D num_hist = singleff_plot->raw_numerator_hist_;
      std::cout << "DEBUG: fitting met bin " << met_bin_idx << std::endl;
      TFitResultPtr frp = eff_plot->data_ratio_plots_[0]->Fit(&linear,"S","",20,50);
      //std::cout << "Fit result, ht bin " << ht_bin_idx << ", met bin " << met_bin_idx << ": " << frp->Parameter(0) << ", " << frp->Parameter(1) << std::endl;
      eff_plot->data_ratio_plots_[0]->Draw("AP");
      linear.Draw("same");
      c->Update();
      c->SaveAs(("plots/trig__sys_0l_fit_htbin"+std::to_string(ht_bin_idx)+"_metbin"+std::to_string(met_bin_idx)+"_"+options.year_string+".pdf").c_str());
      std::cout << "DEBUG: getting fit parameter" << std::endl;
      //double extrapolated_value = frp->Parameter(1);
      double extrapolated_value = linear.GetParameter(1);
      if (extrapolated_value > 1.0) extrapolated_value = 1.0;
      if (extrapolated_value < 0.0) extrapolated_value = 0.0;
      std::cout << "DEBUG: got fit parameter: " << extrapolated_value << std::endl;
      extrapolated_values[ht_bin_idx].push_back(extrapolated_value);
      //bins in 5 GeV, so 20-35 GeV is bins 5-7
      double nominal_value = num_hist.Integral(1,7)/den_hist.Integral(1,7);
      nominal_values[ht_bin_idx].push_back(nominal_value);
      extrapolation_syst[ht_bin_idx].push_back(TMath::Abs(nominal_value-extrapolated_value)/nominal_value);
    }
  }

  //make latex table and vector to use with efficiencies
  ofstream table_file;
  std::vector<std::vector<double>> variations_table;
  std::vector<double> variations_row;
  std::vector<std::vector<double>> ambb_table;
  std::vector<double> ambb_row;
  std::vector<std::vector<double>> reference_table;
  std::vector<double> reference_row;
  std::vector<std::vector<double>> nel_table;
  std::vector<double> nel_row;
  std::vector<std::vector<double>> datamc_table;
  std::vector<double> datamc_row;
  std::vector<std::vector<double>> greater_variations_rows;
  std::vector<std::vector<double>> systematics_rows;
  table_file.open("tables/trig__sys_extrapolation_"+options.year_string+".tex");
  table_file << "\\begin{tabular}{lccccccc}\\hline \\hline \n";
  table_file << "p$_\\text{T}^\\text{miss}$ range& [150,160]& [160,180]& [180,200]& [200,225]& [225,250]& [250,300]& [300+]\\\\ \\hline \\hline \n";
  std::vector<std::string> ht_bin_names = {"$0 \\leq H_{T} < 300$~GeV","$300 \\leq H_{T} < 400$~GeV","$400 \\leq H_{T} < 600$~GeV","$600 \\leq H_{T} < 800$~GeV","$H_{T} \\geq 800$~GeV"};
  std::ostringstream str_stream;
  str_stream << std::fixed;
  for (unsigned int ht_bin_idx = 0; ht_bin_idx < (sys_ht_bins.size()-1); ht_bin_idx++) {
    table_file << ht_bin_names[ht_bin_idx] << " (Nominal)";
    for (unsigned int met_bin_idx = 0; met_bin_idx < (sys_met_bins.size()-1); met_bin_idx++) {
      str_stream << std::setprecision(2);
      str_stream << nominal_values[ht_bin_idx][met_bin_idx];
      table_file << "& " << str_stream.str();
      str_stream.str("");
      str_stream.clear();
    }
    table_file << "\\\\ \n";
    table_file << ht_bin_names[ht_bin_idx] << " Extrapolated";
    for (unsigned int met_bin_idx = 0; met_bin_idx < (sys_met_bins.size()-1); met_bin_idx++) {
      str_stream << std::setprecision(2);
      str_stream << extrapolated_values[ht_bin_idx][met_bin_idx];
      table_file << "& " << str_stream.str();
      str_stream.str("");
      str_stream.clear();
    }
    table_file << "\\\\ \\hline \n";
  }
  table_file << "\\end{tabular} \n";
  table_file.close();
  std::cout << "Wrote tables/trig_sys_extrapolation_"+options.year_string+".tex" << std::endl;
  delete c;
  return extrapolation_syst;
}

std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name) {
  std::string row_text = name;
  std::ostringstream str_stream;
  str_stream << std::fixed;
  double x, y;
  for (int i = 0; i<7; i++) {
    row_text = row_text + "& ";
    gr->GetPoint(i,x,y);
    y = y*100;
    str_stream << std::setprecision(1);
    str_stream << (y);
    row_text = row_text + str_stream.str();
    str_stream.str("");
    str_stream.clear();
    str_stream << std::setprecision(1);
    str_stream << (gr->GetErrorYlow(i)*100);
    row_text = row_text + "$_{-" + str_stream.str();
    str_stream.str("");
    str_stream.clear();
    str_stream << (gr->GetErrorYhigh(i)*100);
    row_text = row_text + "}^{+" + str_stream.str() + "}$";
    str_stream.str("");
    str_stream.clear();
  }
  row_text = row_text + "\\\\";
  return row_text;
}

std::vector<double> get_systematic_table_values(TGraphAsymmErrors *gr) {
  std::vector<double> systematic_table_values;
  double x, y;
        for (int i = 0; i<7; i++) {
    gr->GetPoint(i,x,y);
    y = y*100;
    systematic_table_values.push_back(y);
  }
  return systematic_table_values;
}

std::vector<double> get_systematic_table_variations(std::vector<std::vector<double>> sys_table) {
  //0 is nominal
  std::vector<double> systematic_table_variations;
  for (unsigned int met_idx = 0; met_idx < sys_table[0].size(); met_idx++) {
    double max_difference = 0;
    for (unsigned int row = 1; row < sys_table.size(); row++) {
      if (TMath::Abs(sys_table[0][met_idx]-sys_table[row][met_idx]) > max_difference) {
        max_difference = TMath::Abs(sys_table[0][met_idx]-sys_table[row][met_idx]);
      }
    }
    systematic_table_variations.push_back(max_difference/sys_table[0][met_idx]);
  }
  return systematic_table_variations;
}

std::string make_systematic_table_row(std::vector<double> values, bool total) {
  std::string row_text = "Syst. uncertainty ";
  if (total) {
    row_text = "{\\bf Total Syst. uncertainty} ";
  }
  std::ostringstream str_stream;
  str_stream << std::fixed;
        for (int i = 0; i<7; i++) {
    row_text = row_text + "& ";
    str_stream << std::setprecision(1);
    str_stream << (values[i]*100);
    row_text = row_text + str_stream.str();
    str_stream.str("");
    str_stream.clear();
  }
  row_text = row_text + "\\\\ \\hline \\hline";
  return row_text;
}

std::vector<double> add_variations_quadrature(std::vector<std::vector<double>> values) {
  std::vector<double> variations_quadrature;
  for (unsigned int met_idx = 0; met_idx < values[0].size(); met_idx++) {
    double quad_sum = 0;
    for (unsigned int row = 0; row < values.size(); row++) {
      quad_sum += values[row][met_idx]*values[row][met_idx];
    }
    variations_quadrature.push_back(TMath::Sqrt(quad_sum));
  }
  return variations_quadrature;
}

std::vector<double> max_variations(std::vector<std::vector<double>> values) {
  std::vector<double> r_max_variations;
  for (unsigned int met_idx = 0; met_idx < values[0].size(); met_idx++) {
    double max_var = 0;
    for (unsigned int row = 0; row < values.size(); row++) {
      if (values[row][met_idx] > max_var) {
        max_var = values[row][met_idx];
      }
    }
    r_max_variations.push_back(max_var);
  }
  return r_max_variations;
}


//---------------------------------------------------------------
// efficiency functions
//---------------------------------------------------------------
void trig_postprocess_efficiencies(PlotMaker* pm, 
                                   script_utilities::ArgStruct options,
                                   std::vector<double> systematics, 
                                   std::vector<std::vector<double>> extrapolation_systematics) {
  ofstream apply_effs_file;
  apply_effs_file.open(("src/higgsino/apply_trigeffs_"+options.year_string+".cpp.txt").c_str());
  apply_effs_file << "#include <vector>\n";
  apply_effs_file << "#include \"core/baby.hpp\"\n";
  apply_effs_file << "#include \"core/process.hpp\"\n";
  apply_effs_file << "#include \"core/named_func.hpp\"\n";
  apply_effs_file << "#include \"higgsino/hig_functions.hpp\"\n";
  apply_effs_file << "#include \"higgsino/hig_utilities.hpp\"\n";
  apply_effs_file << "#include \"higgsino/apply_trigeffs" << options.year_string << ".hpp\"\n\n";
  apply_effs_file << "namespace Higfuncs{\n\n";
  //make efficiencies
  //true met
  EfficiencyVariable search_met;
  search_met.name = "met";
  search_met.expression = "b.met()";
  search_met.variable_bins = search_met_bins;
  EfficiencyVariable search_ht;
  search_ht.name = "ht";
  search_ht.expression = "b.ht()";
  search_ht.bins = search_ht_bins;
  generate_2d_efficiency_variable_binning(pm, "get_0l_trigeff"+options.year_string, 
                                          {"trig__eff_search_htbin","_"+options.year_string}, 
                                          search_met, search_ht, apply_effs_file, 
                                          systematics, extrapolation_systematics);
  //fake met 
  EfficiencyVariable qcd_met;
  qcd_met.name = "met";
  qcd_met.expression = "b.met()";
  qcd_met.variable_bins = qcd_met_bins;
  EfficiencyVariable qcd_ht;
  qcd_ht.name = "ht";
  qcd_ht.expression = "b.ht()";
  qcd_ht.bins = qcd_ht_bins;
  generate_2d_efficiency_variable_binning(pm, "get_0l_fakemet_trigeff"+options.year_string, 
                                          {"trig__eff_qcd_htbin","_"+options.year_string}, 
                                          qcd_met, qcd_ht, apply_effs_file, 
                                          {}, {});
  if (HigUtilities::is_in_string_options(options.string_options,"cr")) {
    //1e CR
    EfficiencyVariable onelep_met;
    onelep_met.name = "met";
    onelep_met.expression = "b.met()";
    onelep_met.bins = singlelep_met_bins;
    EfficiencyVariable onelep_ht;
    onelep_ht.name = "ht";
    onelep_ht.expression = "b.ht()";
    onelep_ht.bins = singlelep_ht_bins;
    EfficiencyVariable oneel_pt;
    oneel_pt.name = "el_pt";
    oneel_pt.expression = "Higfuncs::lead_signal_lepton_pt";
    oneel_pt.bins = el_pt_bins;
    generate_3d_efficiency(pm, "get_1el_trigeff"+options.year_string, 
                            {"trig__eff_1e_htbin","_metbin","_"+options.year_string},
                            oneel_pt, onelep_met, onelep_ht, apply_effs_file);
    //1mu CR
    EfficiencyVariable onemu_pt;
    onemu_pt.name = "mu_pt";
    onemu_pt.expression = "Higfuncs::lead_signal_lepton_pt";
    onemu_pt.bins = mu_pt_bins;
    generate_3d_efficiency(pm, "get_1mu_trigeff"+options.year_string, 
                            {"trig__eff_1mu_htbin","_metbin","_"+options.year_string},
                            onemu_pt, onelep_met, onelep_ht, apply_effs_file);
    //2e CR 
    EfficiencyVariable twoel_pt;
    twoel_pt.name = "el_pt";
    twoel_pt.expression = "Higfuncs::lead_signal_lepton_pt";
    twoel_pt.bins = twoel_pt_bins;
    generate_1d_efficiency(pm, "get_2el_trigeff"+options.year_string, "trig__eff_2e_"+options.year_string,
                            twoel_pt, apply_effs_file);
    //2mu CR
    EfficiencyVariable twomu_pt;
    twomu_pt.name = "mu_pt";
    twomu_pt.expression = "Higfuncs::lead_signal_lepton_pt";
    twomu_pt.bins = twomu_pt_bins;
    generate_1d_efficiency(pm, "get_2mu_trigeff"+options.year_string, "trig__eff_2mu_"+options.year_string,
                            twomu_pt, apply_effs_file);
  }
  apply_effs_file << "}\n";
  apply_effs_file.close();
  std::cout << "Wrote src/higgsino/apply_trigeffs_"+options.year_string+".cpp.txt" << std::endl;
}

void generate_3d_efficiency(PlotMaker* pm, std::string func_name, 
                            std::vector<std::string> plot_name,
                            EfficiencyVariable x_var, EfficiencyVariable y_var, 
                            EfficiencyVariable z_var, ofstream &out_file) {
  if (plot_name.size() < 3) {
    ERROR("Error: 3d efficiency requires a 3d array of EfficiencyPlots.");
    return;
  }
  //loop over ht - each ht bin is a new TGraphAsymmErrors
  std::vector<std::vector<TGraphAsymmErrors*>> raw_eff_plots;
  for (unsigned int z_bin = 0; z_bin < z_var.bins.size()-1; z_bin++) {
    std::vector<TGraphAsymmErrors*> raw_eff_plots_row;
    for (unsigned int y_bin = 0; y_bin < y_var.bins.size()-1; y_bin++) {
      EfficiencyPlot * eff_plot = static_cast<EfficiencyPlot*>(pm->GetFigure(plot_name[0]+std::to_string(z_bin)+plot_name[1]+std::to_string(y_bin)+plot_name[2]).get());
      raw_eff_plots_row.push_back(static_cast<TGraphAsymmErrors*>((eff_plot->data_ratio_plots_[0]).get()->Clone()));
    }
    raw_eff_plots.push_back(raw_eff_plots_row);
  }
  //write out function to file
  out_file << "const NamedFunc " << func_name << "(\"" << func_name << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.; // Not used, but for reference\n";
  out_file << "  float eff = 1., " << x_var.name << " = " << x_var.expression << ", " << y_var.name << " = " << y_var.expression << ", " << z_var.name << " = " << z_var.expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  for (unsigned int x_bin = 0; x_bin < x_var.bins.size()-1; x_bin++) {
    double x_upper = x_bin==(x_var.bins.size()-2) ? 9999 : x_var.bins[x_bin+1];
    for (unsigned int y_bin = 0; y_bin < y_var.bins.size()-1; y_bin++) {
      double y_upper = y_bin==(y_var.bins.size()-2) ? 9999 : y_var.bins[y_bin+1];
      for (unsigned int z_bin = 0; z_bin < z_var.bins.size()-1; z_bin++) {
        double z_upper = z_bin==(z_var.bins.size()-2) ? 9999 : z_var.bins[z_bin+1];
        double* eff_pts = raw_eff_plots[z_bin][y_bin]->GetY();
        if (first) {
          out_file << "  if (" << z_var.name << "> ";
          first = false;
        }
        else out_file << "  else if (" << z_var.name << "> ";
        out_file << z_var.bins[z_bin] << " && " << z_var.name << "<= " << z_upper;
        out_file << " && " << y_var.name << "> " << y_var.bins[y_bin] << " && " << y_var.name << "<= " << y_upper;
        out_file << " && " << x_var.name << "> " << x_var.bins[x_bin] << " && " << x_var.name << "<= " << x_upper;
        float errup = raw_eff_plots[z_bin][y_bin]->GetErrorYhigh(x_bin);
        float errdown = raw_eff_plots[z_bin][y_bin]->GetErrorYlow(x_bin);
        out_file << ") {eff = " << eff_pts[x_bin] << "; errup = " << errup << "; errdown = " << errdown << ";}\n";
      }
    }
  }
  out_file << "  std::vector<double> ret = {eff, errup, errdown};\n";
  out_file << "  return ret;\n";
  out_file << "});\n\n";
  for (unsigned int z_bin = 0; z_bin < z_var.bins.size()-1; z_bin++) {
    for (unsigned int y_bin = 0; y_bin < y_var.bins.size()-1; y_bin++) {
      delete raw_eff_plots[z_bin][y_bin];
    }
  }
}

//generates C++ code for 2d efficiencies
//assumes bins for y_var are uniform but bins for x_var are not
void generate_2d_efficiency_variable_binning(PlotMaker* pm, 
                                             std::string func_name,
                                             std::vector<std::string> plot_name,
                                             EfficiencyVariable x_var, 
                                             EfficiencyVariable y_var,
                                             ofstream &out_file, 
                                             std::vector<double> systematics, 
                                             std::vector<std::vector<double>> extrapolation_systematics) {
  if (plot_name.size() < 2) {
    ERROR("Error: 2d efficiency requires a 2d array of EfficiencyPlots.");
    return;
  }
  //loop over ht - each ht bin is a new TGraphAsymmErrors
  std::vector<TGraphAsymmErrors*> raw_eff_plots;
  for (unsigned int y_bin = 0; y_bin < y_var.bins.size()-1; y_bin++) {
    EfficiencyPlot * eff_plot = static_cast<EfficiencyPlot*>(pm->GetFigure(plot_name[0]+std::to_string(y_bin)+plot_name[1]).get());
    raw_eff_plots.push_back(static_cast<TGraphAsymmErrors*>((eff_plot->data_ratio_plots_[0]).get()->Clone()));
  }
  //write out function to file
  out_file << "const NamedFunc " << func_name << "(\"" << func_name << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.;\n";
  out_file << "  float eff = 1., " << x_var.name << " = " << x_var.expression << ", " << y_var.name << " = " << y_var.expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  for (unsigned int y_bin = 0; y_bin < y_var.bins.size()-1; y_bin++) {
    double y_upper = y_bin==(y_var.bins.size()-2) ? 9999 : y_var.bins[y_bin+1];
    double* eff_pts = raw_eff_plots[y_bin]->GetY();
    for (unsigned int x_bin = 0; x_bin < x_var.variable_bins[y_bin].size()-1; x_bin++) {
      double x_upper = x_bin==(x_var.variable_bins[y_bin].size()-2) ? 9999 : x_var.variable_bins[y_bin][x_bin+1];
      if (first) {
        out_file << "  if (" << y_var.name << "> ";
        first = false;
      }
      else out_file << "  else if (" << y_var.name << "> ";
      out_file << y_var.bins[y_bin] << " && " << y_var.name << "<= " << y_upper;
      out_file << " && " << x_var.name << "> " << x_var.variable_bins[y_bin][x_bin] << " && " << x_var.name << "<= " << x_upper;
      float errup = raw_eff_plots[y_bin]->GetErrorYhigh(x_bin);
      float errdown = raw_eff_plots[y_bin]->GetErrorYlow(x_bin);
      //figure out systematics bin and add errors in quadrature, making sure not to exceed 1 or 0
      //sys_met_bins and sys_ht_bins currently hard-coded
      if (systematics.size()>0) {
        for (unsigned int sys_bin = 0; sys_bin < (sys_met_bins.size()-1); sys_bin++) {
          if (sys_met_bins[sys_bin] <= x_var.variable_bins[y_bin][x_bin] && x_var.variable_bins[y_bin][x_bin] < sys_met_bins[sys_bin+1]) {
            errup = TMath::Sqrt(errup*errup+eff_pts[x_bin]*eff_pts[x_bin]*systematics[sys_bin]*systematics[sys_bin]);
            errdown = TMath::Sqrt(errdown*errdown+eff_pts[x_bin]*eff_pts[x_bin]*systematics[sys_bin]*systematics[sys_bin]);
          }
        }
        //apply extrapolation systematic
        for (unsigned int sys_ht_bin = 0; sys_ht_bin < (sys_ht_bins.size()-1); sys_ht_bin++) {
          for (unsigned int sys_met_bin = 0; sys_met_bin < (sys_met_bins.size()-1); sys_met_bin++) {
            if (sys_ht_bins[sys_ht_bin] <= y_var.bins[y_bin] && y_var.bins[y_bin] < sys_ht_bins[sys_ht_bin+1]) {
              if (sys_met_bins[sys_met_bin] <= x_var.variable_bins[y_bin][x_bin] && x_var.variable_bins[y_bin][x_bin] < sys_met_bins[sys_met_bin+1]) {
                errup = TMath::Sqrt(errup*errup+eff_pts[x_bin]*eff_pts[x_bin]*extrapolation_systematics[sys_ht_bin][sys_met_bin]*extrapolation_systematics[sys_ht_bin][sys_met_bin]);
                errdown = TMath::Sqrt(errdown*errdown+eff_pts[x_bin]*eff_pts[x_bin]*extrapolation_systematics[sys_ht_bin][sys_met_bin]*extrapolation_systematics[sys_ht_bin][sys_met_bin]);
              }
            }
          }
        }
        if (errup > (1.0-eff_pts[x_bin])) errup = 1.0-eff_pts[x_bin];
        if (errdown > eff_pts[x_bin]) errdown = eff_pts[x_bin];
      }
      out_file << ") {eff = " << eff_pts[x_bin] << "; errup = " << errup << "; errdown = " << errdown << ";}\n";
    }
  }
  out_file << "  std::vector<double> ret = {eff, errup, errdown};\n";
  out_file << "  return ret;\n";
  out_file << "});\n\n";
  for (TGraphAsymmErrors* raw_eff_plot : raw_eff_plots) {
    delete raw_eff_plot;
  }
}

void generate_1d_efficiency(PlotMaker* pm, std::string func_name, std::string plot_name,
                            EfficiencyVariable x_var, ofstream &out_file) {
  EfficiencyPlot * eff_plot = static_cast<EfficiencyPlot*>(pm->GetFigure(plot_name).get());
  TGraphAsymmErrors* raw_eff_plot = static_cast<TGraphAsymmErrors*>((eff_plot->data_ratio_plots_[0]).get()->Clone());
  //write out functions to file
  //func_name should include year in it
  out_file << "const NamedFunc " << func_name << "(\"" << func_name << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.; // Not used, but for reference\n";
  out_file << "  float eff = 1., " << x_var.name << " = " << x_var.expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  double* eff_x_vals = raw_eff_plot->GetY();
  for (unsigned int x_bin = 0; x_bin < x_var.bins.size()-1; x_bin++) {
    double x_upper = x_bin==(x_var.bins.size()-2) ? 9999 : x_var.bins[x_bin+1];
    if (first) {
      out_file << "  if (" << x_var.name << "> ";
      first = false;
    }
    else out_file << "  else if (" << x_var.name << "> ";
    out_file << x_var.bins[x_bin] << " && " << x_var.name << "<= " << x_upper;
    out_file << ") {eff = " << eff_x_vals[x_bin] << "; errup = " << raw_eff_plot->GetErrorYhigh(x_bin) << "; errdown = " << raw_eff_plot->GetErrorYlow(x_bin) << ";}\n";
  }
  out_file << "  std::vector<double> ret = {eff, errup, errdown};\n";
  out_file << "  return ret;\n";
  out_file << "});\n\n";
  delete raw_eff_plot;
}
