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

using namespace std;
using namespace PlotOptTypes;

namespace{
  //argument options
  bool single_thread = false;
  bool do_variables = false;
  bool do_systematics = false;
  bool do_efficiency = false;
  std::string year_string = "2016";
  std::string out_filename = "triggereff";
  bool do_controlregions = false; //false will omit plots of 1l and 2l CRs, speeding up processing by about 10~20x

  //constants
  const std::vector<double> true_met_bins{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350};
  const std::vector<double> ht_bins{0,250,350,450,600,800,1000,1200};
  const std::vector<double> sys_met_bins{150,160,180,200,225,250,300,350};
  //const std::vector<double> ht_bins{0,200,600,800,1000,1200};
  const std::vector<double> fake_met_bins{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,350,400,450,500,600};
  const std::vector<double> fake_met_bins_ht0to350{150,160,170,180,190,200,225,250,600};
  const std::vector<double> fake_met_bins_ht350to450{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,300,600};
  const std::vector<double> fake_met_bins_ht450to550{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,300,400,600};
  const std::vector<double> fake_met_bins_ht550to650{150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,275,300,400,600};
  const std::vector<double> ht_true_met_bins{0,200,300,400,600,950,1200};
  const std::vector<double> ht_fake_met_bins{0,350,450,550,650,800,1000,1200};
  //const std::vector<double> twodim_met_bins{0,110,120,130,140,150,160,170,180,190,200,210,250};
  //const std::vector<double> el_pt_bins{20,25,30,110,120,150};
  //const std::vector<double> mu_pt_bins{20,25,30,50,100};
  const std::vector<double> twodim_met_bins{0,50,75,100,125,150,175,215,300};
  const std::vector<double> singlelep_ht_bins{0,400,600,800};
  const std::vector<double> el_pt_bins{20,25,30,40,110,120,150};
  const std::vector<double> mu_pt_bins{20,25,30,50,100};
  const std::vector<double> twoel_pt_bins{40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,120};
  const std::vector<double> twomu_pt_bins{40,45,50,60};
  const std::vector<double> true_met_bins_ht0to200 = {150, 155, 160, 165, 170, 180, 190, 200, 350};
  const std::vector<double> true_met_bins_ht200to300 = {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 350};
  const std::vector<double> true_met_bins_ht300to400 = {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 350};
  const std::vector<double> true_met_bins_ht400to600 = {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 350};
  const std::vector<double> true_met_bins_ht600to950 = {150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 350};
  const std::vector<double> true_met_bins_ht950toinf = {150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 275, 300, 350};
  const std::vector<double> sys_ht_bins = {0, 300, 400, 600, 800, 9999};
}

//plot helper function declarations
void generate_ratio_plots(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, const char* hist_title, string hist_name);
void generate_2d_efficiencies(PlotMaker* pm, int denominator_index, int numerator_index, const char* hist_title, string hist_name);
void make_ambb_plots(PlotMaker* pm, int den_index_0, int num_index_0, int den_index_100, int num_index_100, int den_index150, int num_index150, shared_ptr<Process> proc, std::string lumi_string, std::string description, std::string name);
void make_efficiency_plot_wrapper(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, std::string lumi_string, std::string title, std::string xaxis, std::string yaxis, std::string units, std::string out_filename, bool var_width_bins=false, bool is_simulation=false);
void make_efficiency_plot(TH1D* hist_num, TH1D* hist_den, TGraphAsymmErrors* hist_ratio, std::string lumi_string, std::string title, std::string xaxis, std::string yaxis, std::string units, std::string out_filename, bool var_width_bins=false, bool is_simulation=false);

//systematic helper function declarations
std::vector<double> trig_postprocess_systematics(int year, PlotMaker* pm, int first_index, shared_ptr<Process> proc);
std::vector<std::vector<double>> generate_syst_extrapolation_table(int year, PlotMaker* pm, int first_index, shared_ptr<Process> proc);
std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name);
std::vector<double> get_systematic_table_values(TGraphAsymmErrors *gr);
std::vector<double> get_systematic_table_variations(std::vector<std::vector<double>> sys_table);
std::string make_systematic_table_row(std::vector<double> values, bool total=false);
std::vector<double> add_variations_quadrature(std::vector<std::vector<double>> values);
std::vector<double> max_variations(std::vector<std::vector<double>> values);

//efficiency helper function declarations
void trig_postprocess_efficiencies(PlotMaker* pm, int first_index, 
    shared_ptr<Process> met_proc, shared_ptr<Process> lep_proc, std::string year, 
    std::string img_file_extension, std::vector<double> systematics, 
    std::vector<std::vector<double>> extrapolation_systematics);
void draw_3d_efficiency_graphs(PlotMaker* pm, std::vector<int> den_idx, 
    std::vector<int> num_idx, std::vector<double> z_bins, 
    std::vector<double> y_bins, std::vector<double> x_bins, 
    std::string z_var_name, std::string z_title_name, 
    std::string z_expression, std::string y_var_name, std::string y_title_name, 
    std::string y_expression, std::string x_var_name, std::string x_title_name, 
    std::string x_expression, std::string cuts_string, std::string lumi_string, 
    std::string new_trigger_string, std::string plot_name, ofstream &out_file, 
    std::string func_name, std::string year, std::string img_file_extension);
void draw_2d_efficiency_graphs(PlotMaker* pm, int den_idx, int num_idx, 
    std::vector<double> y_bins, std::vector<double> x_bins, std::string y_var_name, 
    std::string y_title_name, std::string y_expression, std::string x_var_name, 
    std::string x_title_name, std::string x_expression, std::string cuts_string, 
    std::string lumi_string, std::string new_trigger_string, std::string plot_name, 
    ofstream &out_file, std::string func_name, std::string year, 
    std::string img_file_extension, std::vector<double> systematics);
void draw_2d_efficiency_graphs_variable_binning(PlotMaker* pm, 
    std::vector<int> den_idx, std::vector<int> num_idx, shared_ptr<Process> proc, 
    std::vector<double> y_bins, std::string y_var_name, std::string y_title_name, 
    std::string y_expression, std::string x_var_name, std::string x_title_name, 
    std::string x_expression, std::string cuts_string, std::string lumi_string, 
    std::string new_trigger_string, std::string plot_name, ofstream &out_file, 
    std::string func_name, std::string year, std::string img_file_extension,
    std::vector<double> systematics, 
    std::vector<std::vector<double>> extrapolation_systematics);
void draw_1d_efficiency_graphs(PlotMaker* pm, int den_idx, int num_idx, 
    shared_ptr<Process> proc, std::vector<double> x_bins, std::string x_var_name, 
    std::string x_title_name, std::string x_expression, std::string cuts_string, 
    std::string lumi_string, std::string new_trigger_string, std::string plot_name, 
    ofstream &out_file, std::string func_name, std::string year, 
    std::string img_file_extension);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  // Define 1D+2D plot types of interest
  Palette colors("txt/colors.txt", "default");
  PlotOpt lin_lumi("txt/plot_styles.txt", "Std1D");
  lin_lumi.Title(TitleType::info).Overflow(OverflowType::overflow); //include overflow to get efficiencies right
  vector<PlotOpt> all_plot_types = {lin_lumi};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).Overflow(OverflowType::overflow)};

  //------------------------------------------------------------------------------------
  //                                 samples
  //------------------------------------------------------------------------------------

  //year-based settings, default to 2016 but changed below
  double lumi = 35.9;
  int year = 2016;
  string data_dir = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2016/data/";
  string mc_dir = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2016/mc/";
  string singleelectron_name = "SingleElectron";

  if (year_string=="2016") year = 2016;
  else if (year_string=="2017") year = 2017;
  else year = 2018;

  vector<shared_ptr<Process> > procs_met150 = {};
  vector<shared_ptr<Process> > procs_met150_mc = {};
  vector<shared_ptr<Process> > procs_1l2j = {};
  shared_ptr<Process> pro_met150;
  shared_ptr<Process> pro_met150_mc;
  shared_ptr<Process> pro_1l2j;

  if (year == 2017) {
  	lumi = 41.5;
  	data_dir = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2017/data/";
  	mc_dir = "/net/cms25/cms25r0/pico/NanoAODv5/higgsino_humboldt/2017/mc/";
  }
  else if (year == 2018) {
  	lumi = 60.0;
  	data_dir = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2018/data/";
  	mc_dir = "/net/cms25/cms25r0/pico/NanoAODv5/higgsino_humboldt/2018/mc/";
	singleelectron_name = "EGamma";
  }
  set<string> str_met150({data_dir+"skim_met150/raw_pico_met150_"+singleelectron_name+"*.root",data_dir+"skim_met150/raw_pico_met150_MET*.root",data_dir+"skim_met150/raw_pico_met150_SingleMuon*.root",data_dir+"skim_met150/raw_pico_met150_JetHT*.root"});
  //set<string> str_met150({data_dir+"skim_met150/raw_pico_met150_*Run2016C*.root"});
  pro_met150 = Process::MakeShared<Baby_pico>("Data MET150", Process::Type::data, kBlack, str_met150, "stitch");
  set<string> str_met150_mc({mc_dir+"skim_met150/pico_met150_TTJets_SingleLeptFromT*.root"});
  pro_met150_mc = Process::MakeShared<Baby_pico>("TTbar 1l", Process::Type::background, kBlack, str_met150_mc, "stitch");
  set<string> str_1l2j({data_dir+"skim_1l2j/raw_pico_1l2j_"+singleelectron_name+"*.root",data_dir+"skim_1l2j/raw_pico_1l2j_MET*.root",data_dir+"skim_1l2j/raw_pico_1l2j_SingleMuon*.root",data_dir+"skim_1l2j/raw_pico_1l2j_JetHT*.root"});
  pro_1l2j = Process::MakeShared<Baby_pico>("Data", Process::Type::data, kBlack, str_1l2j, "stitch");
  procs_met150.push_back(pro_met150);
  procs_met150_mc.push_back(pro_met150_mc);
  procs_1l2j.push_back(pro_1l2j);

  //new version of sample loading. Eventually this should completely replace 
  //the above code
  string base_folder = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/";
  string data_skim_met150_folder = "data/skim_met150/";
  string data_skim_1l2j_folder = "data/skim_1l2j/";
  string mc_skim_met150_folder = "mc/skim_met150/";

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

  std::set<int> years;
  HigUtilities::parseYears(year_string, years);
  std::string total_luminosity_string = HigUtilities::getLuminosityString(year_string);
  std::vector<std::shared_ptr<Process>> procs_data_met150_allyears = {};
  std::vector<std::shared_ptr<Process>> procs_mc_met150_allyears = {};
  std::vector<std::shared_ptr<Process>> procs_data_1l2j_allyears = {};
  std::vector<std::vector<std::shared_ptr<Process>>> procs_data_ambb = {};
  std::vector<std::string> year_strings;
  std::vector<std::string> luminosity_year_string;

  unsigned int year_idx = 0;
  lumi = 0.0;
  for (int iyear : years) {
    std::string iyear_string = "2018";
    std::string iyear_lumi = "60";
    std::set<int> iyear_set = {iyear};
    short iyear_color = kRed;
    if (iyear==2016) { 
      iyear_string = "2016";
      iyear_lumi = "35.9";
      iyear_color = kBlue;
      lumi += 35.9;
    }
    else if (iyear==2017) {
      iyear_string = "2017";
      iyear_lumi = "41.5";
      iyear_color = kGreen;
      lumi += 41.5;
    }
    else {
      lumi += 60.0;
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
  if (do_variables) {
	  //-------------6.2 table 1 plots (0l Real MET: MET, HT(low MET), HT(high MET))-------------
  	pm.Push<EfficiencyPlot>(Axis(100, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
  	    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && 
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_met").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  	    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&150<met&&met<=200" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_ht_lowmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  	    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met&&met<=300" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_ht_highmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
	  //-------------6.3 table 1 plots (0l Fake MET: MET, HT(low MET), HT(high MET))-------------
  	pm.Push<EfficiencyPlot>(Axis(100, 150., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
  	    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && Higfuncs::jetht_trigger,
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakemet_met").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  	    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&150<met&&met<=200" && Higfuncs::jetht_trigger,
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakmet_ht_lowmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  	    Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&200<met&&met<=300" && Higfuncs::jetht_trigger,
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lfakmet_ht_highmet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
	  //-------------9.3 table 1 plots (0l Real MET: Nj Nb DRmax, Ambb(not in AN))-------------
  	pm.Push<EfficiencyPlot>(Axis(5, -0.5, 4.5, nb_higgsino, "N_{b}", {}),
  	    Higfuncs::pass_filters && "njet>=4&&njet<=5&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_nb").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(5, 2.5, 7.5, "njet", "N_{j}", {}),
  	    Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_njet").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(20, 0.0, 4.0, "hig_cand_drmax[0]", "#Delta R_{max}", {}),
  	    Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_drmax").YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
  	pm.Push<EfficiencyPlot>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m_{bb} #GT [GeV]", {}),
  	    Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&200<met" &&
        (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
        Higfuncs::met_trigger,
        procs_data_met150_allyears).Tag("FixName:trig_eff_0lrealmet_ambb");
	  //-------------9.3 table 2 plots (0l Real MET: <mbb> dependence)-------------
    for (unsigned int year_index = 0; year_index < procs_data_ambb.size(); year_index++) {
  	  pm.Push<EfficiencyPlot>(Axis(30, 150., 300., "met", "Offline p_{T}^{miss} [GeV]", {}),
  	      Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300" &&
          (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
          Higfuncs::met_trigger,
          procs_data_ambb[year_index]).Tag("FixName:trig_eff_0lrealmet_met_ht250to300_ambbcuts_"+year_strings[year_index]).YTitle("p_{T}^{miss} Triggers").LuminosityTag(luminosity_year_string[year_index]);
  	  pm.Push<EfficiencyPlot>(Axis(30, 150., 300., "met", "Offline p_{T}^{miss} [GeV]", {}),
  	      Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400" &&
          (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
          Higfuncs::met_trigger,
          procs_data_ambb[year_index]).Tag("FixName:trig_eff_0lrealmet_met_ht300to400_ambbcuts_"+year_strings[year_index]).YTitle("p_{T}^{miss} Triggers").LuminosityTag(luminosity_year_string[year_index]);
  	  pm.Push<EfficiencyPlot>(Axis(30, 150., 300., "met", "Offline p_{T}^{miss} [GeV]", {}),
  	      Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600" &&
          (Higfuncs::lead_signal_lepton_pt < 35) && Higfuncs::el_trigger, 
          Higfuncs::met_trigger,
          procs_data_ambb[year_index]).Tag("FixName:trig_eff_0lrealmet_met_ht400to600_ambbcuts_"+year_strings[year_index]).YTitle("p_{T}^{miss} Triggers").LuminosityTag(luminosity_year_string[year_index]);
    }
	  if (do_controlregions) {
		  //-------------6.4 table 1 plots (1l CR: MET lep_pt HT(lowmet) HT(highmet))-------------
      //electrons
  		pm.Push<EfficiencyPlot>(Axis(55, 0., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nel==1" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1el_met").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Electron p_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nel==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1el_leppt").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nel==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1el_ht_lowmet").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nel==1&&200<=met&&met<300" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::el_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1el_ht_highmet").YTitle("p_{T}^{miss} and Electron Triggers").LuminosityTag(total_luminosity_string);
      //muons
  		pm.Push<EfficiencyPlot>(Axis(55, 0., 550., "met", "Offline p_{T}^{miss} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nmu==1" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1mu_met").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Muon p_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nmu==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1mu_leppt").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nmu==1&&150<=met&&met<200" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1mu_ht_lowmet").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(13, 0., 1300., "ht", "H_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nmu==1&&200<=met&&met<300" && Higfuncs::jetht_trigger, 
          (Higfuncs::met_trigger||Higfuncs::mu_trigger),
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_1mu_ht_highmet").YTitle("p_{T}^{miss} and Muon Triggers").LuminosityTag(total_luminosity_string);
		  //-------------6.5 table 1 plots (2l CR: lep_pt HT)-------------
      //electrons
  		pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::el_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_2el_leppt").YTitle("Electron Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::el_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_2el_ht").YTitle("Electron Triggers").LuminosityTag(total_luminosity_string);
      //muons
  		pm.Push<EfficiencyPlot>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::mu_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_2mu_leppt").YTitle("Muon Triggers").LuminosityTag(total_luminosity_string);
  		pm.Push<EfficiencyPlot>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
  		    Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), 
          Higfuncs::mu_trigger,
          procs_data_1l2j_allyears).Tag("FixName:trig_eff_2mu_ht").YTitle("Muon Triggers").LuminosityTag(total_luminosity_string);
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
  if (do_systematics) {
	//nj3 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj3_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj3_sys_num");
	//nj4
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet==4&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj4_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet==4&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj4_sys_num");
	//nj5
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet==5&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj5_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet==5&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj5_sys_num");
	//nb0
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && nb_higgsino==0., procs_met150, all_plot_types).Tag("FixName:trig_raw_nb0_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && nb_higgsino==0. && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb0_sys_num");
	//nb1
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && nb_higgsino==1., procs_met150, all_plot_types).Tag("FixName:trig_raw_nb1_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && nb_higgsino==1. && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb1_sys_num");
	//nb2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && nb_higgsino>=2., procs_met150, all_plot_types).Tag("FixName:trig_raw_nb2_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && nb_higgsino>=2. && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb2_sys_num");
	//drmax < 2.2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]<2.2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax0_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]<2.2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax0_sys_num");
	//drmax > 2.2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]>=2.2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax22_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]>=2.2" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax22_sys_num");
	//300<=ht<400 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_ht300to400_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_ht300to400_sys_num");
	//am<140, 300<=ht<400
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]<140" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am0_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]<140" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am0_sys_num");
	//am>140, 300<=ht<400
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]>=140" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am140_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]>=140" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am140_sys_num");
	//jet_pt>500 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_refel_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && high_pt_jet && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_refel_sys_num");
	//jet_pt>500, ref trigger jet500
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_refjet_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                Higfuncs::pass_filters && "HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && high_pt_jet && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_refjet_sys_num");
	////jet_pt>500 nel=1, ref trigger jet500 (nominal)
  //	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  //	                Higfuncs::pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_1e_sys_den");
  //	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  //	                Higfuncs::pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet && metnoHigfuncs::mu_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_1e_sys_num");
	////jet_pt>500 nmu=1, ref trigger jet500
  //	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  //	                Higfuncs::pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nmu==1" && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_1mu_sys_den");
  //	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  //	                Higfuncs::pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nmu==1" && high_pt_jet && metnoHigfuncs::mu_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_1mu_sys_num");
    //const std::vector<double> sys_met_bins{150,160,180,200,225,250,300,350};
    //const std::vector<double> new_met_bins = {150
    //std::vector<std::vector<double>> new_met_bins = {};
    //new_met_bins.push_back(std::vector<double>({150, 155, 160, 165, 170, 180, 190, 200, 9999}));
    //new_met_bins.push_back(std::vector<double>({150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 9999}));
    //new_met_bins.push_back(std::vector<double>({150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 9999}));
    //new_met_bins.push_back(std::vector<double>({150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 9999}));
    //new_met_bins.push_back(std::vector<double>({150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 275, 300, 9999}));
    //new_met_bins.push_back(std::vector<double>({150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 275, 300, 9999}));
    for (unsigned int ht_bin_idx = 0; ht_bin_idx < (sys_ht_bins.size()-1); ht_bin_idx++) {
      for (unsigned int met_bin_idx = 0; met_bin_idx < (sys_met_bins.size()-1); met_bin_idx++) {
    	  pm.Push<EfficiencyPlot>(Axis(16, 0., 80., Higfuncs::lead_signal_lepton_pt, "Electron p_{T} [GeV]", {}),
    	      Higfuncs::pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && Higfuncs::el_trigger && 
            (ht >= sys_ht_bins[ht_bin_idx]) && (ht < sys_ht_bins[ht_bin_idx+1]) && 
            (met >= sys_met_bins[met_bin_idx]) && (met < sys_met_bins[met_bin_idx+1]),
            Higfuncs::met_trigger,
            procs_met150).Tag("FixName:trig_eff_0l_htbin"+std::to_string(ht_bin_idx)+"_metbin"+std::to_string(met_bin_idx)+"_"+year_string).YTitle("p_{T}^{miss} Triggers").LuminosityTag(total_luminosity_string);
        //TODO: fix to allow generating efficiencies for multiple years at once?
      }
    }
  }

  // trigger efficiency plots and application file
  if (do_efficiency) {
    //0l Real MET HT bins
  	pm.Push<Hist1D>(Axis(true_met_bins_ht0to200, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "0<=ht&&ht<200", 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht0to200_den");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht0to200, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "0<=ht&&ht<200"
                    && Higfuncs::met_trigger, 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht0to200_num");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht200to300, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "200<=ht&&ht<300", 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht200to300_den");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht200to300, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "200<=ht&&ht<300"
                    && Higfuncs::met_trigger, 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht200to300_num");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht300to400, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "300<=ht&&ht<400", 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht300to400_den");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht300to400, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "300<=ht&&ht<400"
                    && Higfuncs::met_trigger, 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht300to400_num");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht400to600, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "400<=ht&&ht<600", 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht400to600_den");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht400to600, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "400<=ht&&ht<600"
                    && Higfuncs::met_trigger, 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht400to600_num");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht600to950, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "600<=ht&&ht<950", 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht600to950_den");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht600to950, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "600<=ht&&ht<950"
                    && Higfuncs::met_trigger, 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht600to950_num");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht950toinf, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "950<=ht", 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht950toInf_den");
  	pm.Push<Hist1D>(Axis(true_met_bins_ht950toinf, "met", "MET", {}),
  	                Higfuncs::pass_filters && "!low_dphi_met&&njet>=3&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && Higfuncs::el_trigger && "950<=ht"
                    && Higfuncs::met_trigger, 
                    procs_met150, all_plot_types).Tag("FixName:trig_raw_0l_eff_ht950toInf_num");
	  //QCD HT bins
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht0to350, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&ht<350" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht0to350_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht0to350, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&ht<350" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht0to350_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht350to450, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&350<=ht&&ht<450" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht350to450_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht350to450, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&350<=ht&&ht<450" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht350to450_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht450to550, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&450<=ht&&ht<550" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht450to550_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht450to550, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&450<=ht&&ht<550" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht450to550_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht550to650, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&550<=ht&&ht<650" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht550to650_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht550to650, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&550<=ht&&ht<650" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht550to650_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&650<=ht&&ht<800" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht650to800_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&650<=ht&&ht<800" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht650to800_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&800<=ht&&ht<1000" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht800to1000_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&800<=ht&&ht<1000" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht800to1000_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&1000<=ht" && Higfuncs::jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht1000toInf_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0&&1000<=ht" && Higfuncs::jetht_trigger && Higfuncs::met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht1000toInf_num");
  	//pm.Push<Hist2D>(Axis(fake_met_bins, "met", "MET", {}),
  	//                Axis(ht_fake_met_bins, "ht", "HT", {}),
  	//                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && ht_trigger, procs_met150, twodim_plotopts).Tag("FixName:trig_raw_qcd_eff_den");
  	//pm.Push<Hist2D>(Axis(fake_met_bins, "met", "MET", {}),
  	//                Axis(ht_fake_met_bins, "ht", "HT", {}),
  	//                Higfuncs::pass_filters && "low_dphi_met&&nvlep==0" && ht_trigger && Higfuncs::met_trigger, procs_met150, twodim_plotopts).Tag("FixName:trig_raw_qcd_eff_num");
	if (do_controlregions) {
		//1e
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==1&&0<=ht&&ht<400" && Higfuncs::jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht0to400_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==1&&0<=ht&&ht<400" && Higfuncs::jetht_trigger && (Higfuncs::met_trigger||Higfuncs::el_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht0to400_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==1&&400<=ht&&ht<600" && Higfuncs::jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht400to600_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==1&&400<=ht&&ht<600" && Higfuncs::jetht_trigger && (Higfuncs::met_trigger||Higfuncs::el_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht400to600_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==1&&600<=ht" && Higfuncs::jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht600toInf_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==1&&600<=ht" && Higfuncs::jetht_trigger && (Higfuncs::met_trigger||Higfuncs::el_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht600toInf_eff_num");
		//1mu
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==1&&0<=ht&&ht<=400" && Higfuncs::jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht0to400_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==1&&0<=ht&&ht<=400" && Higfuncs::jetht_trigger && (Higfuncs::met_trigger||Higfuncs::mu_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht0to400_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==1&&400<=ht&&ht<=600" && Higfuncs::jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht400to600_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==1&&400<=ht&&ht<=600" && Higfuncs::jetht_trigger && (Higfuncs::met_trigger||Higfuncs::mu_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht400to600_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==1&&600<=ht" && Higfuncs::jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht600toInf_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==1&&600<=ht" && Higfuncs::jetht_trigger && (Higfuncs::met_trigger||Higfuncs::mu_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht600toInf_eff_num");
		//2l
  		pm.Push<Hist1D>(Axis(twoel_pt_bins, Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2el_eff_den");
  		pm.Push<Hist1D>(Axis(twoel_pt_bins, Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger) && Higfuncs::el_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2el_eff_num");
  		pm.Push<Hist1D>(Axis(twomu_pt_bins, Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2mu_eff_den");
  		pm.Push<Hist1D>(Axis(twomu_pt_bins, Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                Higfuncs::pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (Higfuncs::met_trigger||Higfuncs::jetht_trigger) && Higfuncs::mu_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2mu_eff_num");
	}
  }

  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.print_2d_figures_ = false;
  pm.MakePlots(lumi);

  //-------------------------------------------------------------------
  // post processing
  //-------------------------------------------------------------------
  std::string str_lumi = "35.9 fb^{-1} (13 TeV)";
  std::string str_mc = "Summer16 t#bar{t} 1l";
  std::string str_met_trig = "MET Triggers";
  std::string str_el_trig = "SingleElectron Triggers";
  std::string str_mu_trig = "SingleMuon Triggers";
  std::string str_realmet_ref_trig = "#it{SingleElectron}";
  std::string str_fakemet_ref_trig = "#it{JetHT}";
  //std::string str_fakemet_ref_trig = "#it{Jet500}";
  if (year == 2017) {
    str_lumi = "41.5 fb^{-1} (13 TeV)";
    str_mc = "Fall17 t#bar{t} 1l";
  }
  else if (year == 2018) {
    str_lumi = "60 fb^{-1} (13 TeV)";
    str_mc = "Autumn18 t#bar{t} 1l";
  }

  //also save plots to root file in case something fails
  TFile* out_file = TFile::Open(("ntuples/"+out_filename+".root").c_str(),"RECREATE");
  int pm_idx = 0;
  if (do_variables) {
    //pm_idx += (10+3*procs_data_ambb.size());
    if (do_controlregions) {
      pm_idx += 4*3;
    }
    //pm_idx += 2;
    pm_idx += 7;
    //pm_idx += 7;
  }
  std::vector<double> systematics;
  std::vector<std::vector<double>> extrapolation_systematics;
  if (do_systematics) {
    systematics = trig_postprocess_systematics(year, &pm, pm_idx, pro_met150);
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj3");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}= 4 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj4");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}= 5 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj5");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 N_{b}=0 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb0");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 N_{b}=1 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb1");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 N_{b}#geq 1 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb2");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 #Delta R_{max}#leq 2.2 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_drmaxlow");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 #Delta R_{max}>2.2 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_drmaxhigh");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 300 #leq HT < 400 GeV, high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fixht");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 300 #leq HT < 400 GeV, #LT m#RT #leq 140 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_amlow");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 300 #leq HT < 400 GeV, #LT m#RT > 140 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_amhigh");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 p_{Tjet}>500 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_isoeljet500");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 3 p_{Tjet}>500 GeV high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_jetjet500");
    pm_idx += 2;
    //generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 3 p_{Tjet}>500 GeV, high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_1e");
    //pm_idx += 2;
    //generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 3 p_{Tjet}>500 GeV, high #Delta#phi N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_1mu");
    //pm_idx += 2;
    std::cout << "DEBUG: starting extrapolation systematics" << std::endl;
    extrapolation_systematics = generate_syst_extrapolation_table(year, &pm, pm_idx, pro_met150);
    pm_idx += (sys_ht_bins.size()-1)*(sys_met_bins.size()-1);
  }
  if (do_efficiency) {
    std::cout << "DEBUG: starting efficiency" << std::endl;
    trig_postprocess_efficiencies(&pm, pm_idx, pro_met150, pro_1l2j, year_string, "pdf", systematics, extrapolation_systematics);
    std::cout << "DEBUG: finishing efficiency" << std::endl;
    //real MET plots
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 p_{Te}<35 GeV 0<H_{T}<200 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV]", "hist_realmet_ht0to200");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 p_{Te}<35 GeV 200<H_{T}<300 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV]", "hist_realmet_ht200to300");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 p_{Te}<35 GeV 300<H_{T}<400 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV]", "hist_realmet_ht300to400");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 p_{Te}<35 GeV 400<H_{T}<600 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV]", "hist_realmet_ht400to600");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 p_{Te}<35 GeV 600<H_{T}<950 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV]", "hist_realmet_ht600to950");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 p_{Te}<35 GeV 900<H_{T} GeV, 35.9 fb^{-1} (13 TeV); MET [GeV]", "hist_realmet_ht950toInf");
    pm_idx += 2;
    //fake MET plots
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 0<H_{T}<350 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht0to350");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 350<H_{T}<450 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht350to450");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 450<H_{T}<550 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht450to550");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 550<H_{T}<650 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht550to650");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 650<H_{T}<800 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht650to800");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 800<H_{T}<1000 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht800to1000");
    pm_idx += 2;
    generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, H_{T}>1000 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemet_eff_ht1000toInf");
    pm_idx += 2;
    //generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_fakemetht");
    //pm_idx += 2;
    if (do_controlregions) {
      generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||Ele(27_WPTight||35_WPTight||115)] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 0#leq HT < 400 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], Electron pt [GeV]", "hist_metelptht0to400");
      pm_idx += 2;
      generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||Ele(27_WPTight||35_WPTight||115)] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 400#leq HT < 600 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], Electron pt [GeV]", "hist_metelptht400to600");
      pm_idx += 2;
      generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||Ele(27_WPTight||35_WPTight||115)] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 600 GeV #leq HT , 35.9 fb^{-1} (13 TeV); MET [GeV], Electron pt [GeV]", "hist_metelptht600toInf");
      pm_idx += 2;
      generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||IsoMu(24||27)||Mu50] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 0#leq HT < 400 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], Muon pt [GeV]", "hist_metmuptht0to400");
      pm_idx += 2;
      generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||IsoMu(24||27)||Mu50] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 400#leq HT < 600 GeV, 35.9 fb^{-1} (13 TeV); MET [GeV], Muon pt [GeV]", "hist_metmuptht400to600");
      pm_idx += 2;
      generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)||IsoMu(24||27)||Mu50] Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 600 GeV #leq HT, 35.9 fb^{-1} (13 TeV); MET [GeV], Muon pt [GeV]", "hist_metmuptht600toInf");
      pm_idx += 2;
      generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{e}=2 80<m_{ee}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Electron Pt [GeV]; Efficiency [Ele(27_WPTight||35_WPTight||115)]","hist_elel");
      pm_idx += 2;
      generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{#mu}=2 80<m_{#mu#mu}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Muon Pt [GeV]; Efficiency [IsoMu(24||27)||Mu50]","hist_mumu");
      pm_idx += 2;
	  }
  }
  out_file->Close();

  time(&endtime);
  cout<<endl<<"Processing took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

//---------------------------------------------------------------
//plot functions
//---------------------------------------------------------------
void generate_ratio_plots(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, const char* hist_title, string hist_name) {
  Hist1D * den_h1d = static_cast<Hist1D*>(pm->Figures()[denominator_index].get());
  Hist1D::SingleHist1D* den_sh = static_cast<Hist1D::SingleHist1D*>(den_h1d->GetComponent(proc.get()));
  TH1D den_h = den_sh->scaled_hist_;
  Hist1D * num_h1d = static_cast<Hist1D*>(pm->Figures()[numerator_index].get());
  Hist1D::SingleHist1D* num_sh = static_cast<Hist1D::SingleHist1D*>(num_h1d->GetComponent(proc.get()));
  TH1D num_h = num_sh->scaled_hist_;
  TGraphAsymmErrors* ratio_h = new TGraphAsymmErrors(&num_h,&den_h,"cp");
  ratio_h->SetTitle(hist_title);
  ratio_h->GetYaxis()->SetRangeUser(0,1.4);
  ratio_h->SetMarkerStyle(kFullCircle);
  ratio_h->SetMarkerSize(1);
  ratio_h->GetXaxis()->SetTitleSize(0.04);
  ratio_h->GetYaxis()->SetTitleSize(0.04);
  ratio_h->GetXaxis()->SetTitleOffset(1.0);
  ratio_h->GetYaxis()->SetTitleOffset(1.0);
  den_h.Write((hist_name+"_denominator").c_str());
  num_h.Write((hist_name+"_numerator").c_str());
  ratio_h->Write((hist_name+"_ratio").c_str());
}

void generate_2d_efficiencies(PlotMaker* pm, int denominator_index, int numerator_index, const char* hist_title, string hist_name) {
  Hist2D * effhist_den = static_cast<Hist2D*>(pm->Figures()[denominator_index].get());
  TH2D eff_den_h = effhist_den->GetDataHist();
  Hist2D * effhist_num = static_cast<Hist2D*>(pm->Figures()[numerator_index].get());
  TH2D eff_num_h = effhist_num->GetDataHist();
  TH2D *eff_ratio_h = static_cast<TH2D*>(eff_num_h.Clone());
  eff_ratio_h->Divide(&eff_den_h);
  eff_ratio_h->SetTitle(hist_title);
  eff_den_h.Write((hist_name+"_denominator").c_str());
  eff_num_h.Write((hist_name+"_numerator").c_str());
  eff_ratio_h->Write((hist_name+"_ratio").c_str());
}

void make_efficiency_plot_wrapper(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, std::string lumi_string, std::string title, std::string xaxis, std::string yaxis, std::string units, std::string out_filename, bool var_width_bins, bool is_simulation) {
  Hist1D * den_h1d = static_cast<Hist1D*>(pm->Figures()[denominator_index].get());
  Hist1D::SingleHist1D* den_sh = static_cast<Hist1D::SingleHist1D*>(den_h1d->GetComponent(proc.get()));
  TH1D den_h = den_sh->scaled_hist_;
  Hist1D * num_h1d = static_cast<Hist1D*>(pm->Figures()[numerator_index].get());
  Hist1D::SingleHist1D* num_sh = static_cast<Hist1D::SingleHist1D*>(num_h1d->GetComponent(proc.get()));
  TH1D num_h = num_sh->scaled_hist_;
  TGraphAsymmErrors* ratio_h = new TGraphAsymmErrors(&num_h,&den_h,"cp");
  ratio_h->GetYaxis()->SetRangeUser(0,1.4);
  ratio_h->SetMarkerStyle(kFullCircle);
  ratio_h->SetMarkerSize(1);
  ratio_h->GetXaxis()->SetTitleSize(0.04);
  ratio_h->GetYaxis()->SetTitleSize(0.04);
  ratio_h->GetXaxis()->SetTitleOffset(1.0);
  ratio_h->GetYaxis()->SetTitleOffset(1.0);
  make_efficiency_plot(&num_h, &den_h, ratio_h, lumi_string, title, xaxis, yaxis, units, out_filename, var_width_bins, is_simulation);
}

void make_ambb_plots(PlotMaker* pm, int den_index_0, int num_index_0, int den_index_100, int num_index_100, int den_index_150, int num_index_150, shared_ptr<Process> proc, std::string lumi_string, std::string description, std::string name) {
  //get the histograms from the plot maker
  Hist1D * den0_h1d = static_cast<Hist1D*>(pm->Figures()[den_index_0].get());
  Hist1D::SingleHist1D* den0_sh = static_cast<Hist1D::SingleHist1D*>(den0_h1d->GetComponent(proc.get()));
  TH1D den0_h = den0_sh->scaled_hist_;
  Hist1D * num0_h1d = static_cast<Hist1D*>(pm->Figures()[num_index_0].get());
  Hist1D::SingleHist1D* num0_sh = static_cast<Hist1D::SingleHist1D*>(num0_h1d->GetComponent(proc.get()));
  TH1D num0_h = num0_sh->scaled_hist_;
  TGraphAsymmErrors* ratio0_h = new TGraphAsymmErrors(&num0_h,&den0_h,"cp");

  Hist1D * den100_h1d = static_cast<Hist1D*>(pm->Figures()[den_index_100].get());
  Hist1D::SingleHist1D* den100_sh = static_cast<Hist1D::SingleHist1D*>(den100_h1d->GetComponent(proc.get()));
  TH1D den100_h = den100_sh->scaled_hist_;
  Hist1D * num100_h1d = static_cast<Hist1D*>(pm->Figures()[num_index_100].get());
  Hist1D::SingleHist1D* num100_sh = static_cast<Hist1D::SingleHist1D*>(num100_h1d->GetComponent(proc.get()));
  TH1D num100_h = num100_sh->scaled_hist_;
  TGraphAsymmErrors* ratio100_h = new TGraphAsymmErrors(&num100_h,&den100_h,"cp");

  Hist1D * den150_h1d = static_cast<Hist1D*>(pm->Figures()[den_index_150].get());
  Hist1D::SingleHist1D* den150_sh = static_cast<Hist1D::SingleHist1D*>(den150_h1d->GetComponent(proc.get()));
  TH1D den150_h = den150_sh->scaled_hist_;
  Hist1D * num150_h1d = static_cast<Hist1D*>(pm->Figures()[num_index_150].get());
  Hist1D::SingleHist1D* num150_sh = static_cast<Hist1D::SingleHist1D*>(num150_h1d->GetComponent(proc.get()));
  TH1D num150_h = num150_sh->scaled_hist_;
  TGraphAsymmErrors* ratio150_h = new TGraphAsymmErrors(&num150_h,&den150_h,"cp");

  //draw ratio and axis labels
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas(("can_"+out_filename).c_str(),"can",1024,1024);
  can->cd();
  TPad* pad = new TPad(("pad_"+out_filename).c_str(),"pad",0.,0.,1.0,1.0);
  pad->Draw();
  pad->cd();
  pad->SetGrid();
  gPad->SetMargin(0.15,0.15,0.15,0.15);

  ratio0_h->GetXaxis()->SetRangeUser(150., 300.);
  ratio0_h->GetYaxis()->SetTitleOffset(1.4);
  ratio0_h->GetYaxis()->SetRangeUser(0, 1.4);
  ratio0_h->SetTitle("; p_{T}^{miss} [GeV]; Efficiency MET Triggers");
  ratio0_h->SetLineWidth(2);
  ratio0_h->SetMarkerColor(kRed);
  ratio0_h->SetLineColor(kRed);
  ratio0_h->Draw("AP");
  ratio100_h->SetLineWidth(2);
  ratio100_h->SetMarkerColor(kGreen);
  ratio100_h->SetLineColor(kGreen);
  ratio100_h->Draw("P");
  ratio150_h->SetLineWidth(2);
  ratio150_h->SetMarkerColor(kBlue);
  ratio150_h->SetLineColor(kBlue);
  ratio150_h->Draw("P");
  
  //draw overlay text
  TLatex t;
  t.SetTextColor(kBlack);
  t.SetTextSize(0.04);
  t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  t.SetTextAlign(31);
  t.DrawLatexNDC(0.845,0.87,("#font[42]{"+lumi_string+"}").c_str());
  t.SetTextAlign(33);
  t.SetTextSize(0.02);
  t.DrawLatexNDC(0.825,0.83,("#font[42]{Denom: #it{SingleElectron}, N_{j}#geq 4, hi-#Delta #phi, N_{e}=1, " +description+"}").c_str());

  //draw legend
  TLegend * legend = new TLegend(0.5, 0.2, 0.8, 0.4);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(ratio0_h,"#LT m #GT < 100 GeV","l");
  legend->AddEntry(ratio100_h,"100 #leq #LT m #GT < 150 GeV","l");
  legend->AddEntry(ratio150_h,"#LT m #GT #geq 150 GeV","l");
  legend->Draw();
  
  //save
  can->Update();
  can->SaveAs(("plots/trig_eff_resolved_metambb_"+name+".pdf").c_str());
}

void make_efficiency_plot(TH1D* hist_num, TH1D* hist_den, TGraphAsymmErrors* hist_ratio, std::string lumi_string, std::string title, std::string xaxis, std::string yaxis, std::string units, std::string out_filename, bool var_width_bins, bool is_simulation) {
	//draw ratio and axis labels
	gStyle->SetOptStat(0);
	TCanvas* can = new TCanvas(("can_"+out_filename).c_str(),"can",1024,1024);
	can->cd();
	TPad* pad = new TPad(("pad_"+out_filename).c_str(),"pad",0.,0.,1.0,1.0);
	pad->Draw();
	pad->cd();
	pad->SetGrid();
	gPad->SetMargin(0.15,0.15,0.15,0.15);

	double xlow = hist_den->GetBinLowEdge(1);
	double xhigh = hist_den->GetBinLowEdge(hist_den->GetNbinsX())+hist_den->GetBinWidth(hist_den->GetNbinsX());
	double bin_size = (xhigh-xlow)/hist_den->GetNbinsX();
	hist_ratio->GetXaxis()->SetRangeUser(xlow, xhigh);
	if (yaxis.size() > 40) {
		hist_ratio->GetYaxis()->SetTitleSize(0.025);
		hist_ratio->GetYaxis()->SetTitleOffset(1.6);
	}
	else {
		hist_ratio->GetYaxis()->SetTitleOffset(1.4);
	}
	hist_ratio->GetYaxis()->SetRangeUser(0, 1.4);
	if (units == "") {
		hist_ratio->SetTitle((";"+xaxis+";"+yaxis).c_str());
	}
	else {
		hist_ratio->SetTitle((";"+xaxis+"("+units+");"+yaxis).c_str());
	}
	hist_ratio->SetLineWidth(2);
	hist_ratio->Draw("AP");

	//draw overlay text
	TLatex t;
	t.SetTextColor(kBlack);
	t.SetTextSize(0.04);
	if (!is_simulation) {
		t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
	}
	else {
		t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
	}
	t.SetTextAlign(31);
	t.DrawLatexNDC(0.845,0.87,("#font[42]{"+lumi_string+"}").c_str());
	t.SetTextAlign(33);
	t.SetTextSize(0.03);
	if (title.size() < 66) {
		t.DrawLatexNDC(0.825,0.83,("#font[42]{"+title+"}").c_str());
	}
	else {
		//title too long, draw on separate lines
		unsigned int right_split_comma_pos = 999;
		for (unsigned int string_pos = 65; string_pos > 0; string_pos --) {
			if (title[string_pos]==',') {
				//split title at this comma
				right_split_comma_pos = string_pos;
				break;
			}
		}
		if (right_split_comma_pos == 999) {
			//unable to find splitting comma, just let it draw off histogram
			t.DrawLatexNDC(0.825,0.83,("#font[42]{"+title+"}").c_str());
		}
		else {
			t.DrawLatexNDC(0.825,0.83,("#font[42]{"+title.substr(0,right_split_comma_pos+1)+"}").c_str());
			t.DrawLatexNDC(0.825,0.78,("#font[42]{"+title.substr(right_split_comma_pos+1,title.size()-right_split_comma_pos-1)+"}").c_str());
		}
	}
	
	//draw numerator and denominator histograms
	double hist_den_max = hist_den->GetBinContent(hist_den->GetMaximumBin());
	hist_den->Scale(0.5/hist_den_max);
	hist_den->SetLineStyle(2);
	hist_den->SetLineColor(kBlue);
	hist_den->SetLineWidth(2);
	hist_den->SetFillStyle(0);
	hist_den->Draw("same hist");

	hist_num->Scale(0.5/hist_den_max);
	hist_num->SetLineColor(kBlue);
	hist_num->SetLineWidth(2);
	hist_num->SetFillStyle(0);
	//hist_num->SetFillColor(kBlue);
	hist_num->Draw("same hist");

	hist_ratio->Draw("P");

	//draw right axis
	pad->Update();
	TGaxis *norm_axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, 2.8*hist_den_max, 505, "+L");
	norm_axis->SetTickLength(0.3);
	norm_axis->SetLabelSize(0.03);
	if (units != "" && !var_width_bins) {
		norm_axis->SetTitle(("Events/("+std::to_string(static_cast<int>(bin_size))+" "+units+")").c_str());
	}
	else {
		norm_axis->SetTitle("Events/bin");
	}
	norm_axis->SetTitleColor(kBlue);
	norm_axis->SetTitleFont(42);
	norm_axis->SetTitleOffset(1.6);
	norm_axis->SetLineColor(kBlue);
	norm_axis->SetLabelColor(kBlue);
	norm_axis->Draw();

	//save
	can->Update();
	can->SaveAs(("plots/trig_eff_resolved_"+out_filename).c_str());
}


//---------------------------------------------------------------
// systematics functions
//---------------------------------------------------------------
std::vector<double> trig_postprocess_systematics(int year, PlotMaker* pm, int first_index, shared_ptr<Process> proc) {
  //generate ratio plots
  //nj3
  Hist1D * den_h_nj3 = static_cast<Hist1D*>(pm->Figures()[first_index].get());
  Hist1D::SingleHist1D* den_sh_nj3 = static_cast<Hist1D::SingleHist1D*>(den_h_nj3->GetComponent(proc.get()));
  TH1D den_nj3 = den_sh_nj3->scaled_hist_;
  Hist1D * num_h_nj3 = static_cast<Hist1D*>(pm->Figures()[first_index+1].get());
  Hist1D::SingleHist1D* num_sh_nj3 = static_cast<Hist1D::SingleHist1D*>(num_h_nj3->GetComponent(proc.get()));
  TH1D num_nj3 = num_sh_nj3->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_nj3 = new TGraphAsymmErrors(&num_nj3,&den_nj3,"cp");
  //nj4
  Hist1D * den_h_nj4 = static_cast<Hist1D*>(pm->Figures()[first_index+2].get());
  Hist1D::SingleHist1D* den_sh_nj4 = static_cast<Hist1D::SingleHist1D*>(den_h_nj4->GetComponent(proc.get()));
  TH1D den_nj4 = den_sh_nj4->scaled_hist_;
  Hist1D * num_h_nj4 = static_cast<Hist1D*>(pm->Figures()[first_index+3].get());
  Hist1D::SingleHist1D* num_sh_nj4 = static_cast<Hist1D::SingleHist1D*>(num_h_nj4->GetComponent(proc.get()));
  TH1D num_nj4 = num_sh_nj4->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_nj4 = new TGraphAsymmErrors(&num_nj4,&den_nj4,"cp");
  //nj5
  Hist1D * den_h_nj5 = static_cast<Hist1D*>(pm->Figures()[first_index+4].get());
  Hist1D::SingleHist1D* den_sh_nj5 = static_cast<Hist1D::SingleHist1D*>(den_h_nj5->GetComponent(proc.get()));
  TH1D den_nj5 = den_sh_nj5->scaled_hist_;
  Hist1D * num_h_nj5 = static_cast<Hist1D*>(pm->Figures()[first_index+5].get());
  Hist1D::SingleHist1D* num_sh_nj5 = static_cast<Hist1D::SingleHist1D*>(num_h_nj5->GetComponent(proc.get()));
  TH1D num_nj5 = num_sh_nj5->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_nj5 = new TGraphAsymmErrors(&num_nj5,&den_nj5,"cp");
  //nb0
  Hist1D * den_h_nb0 = static_cast<Hist1D*>(pm->Figures()[first_index+6].get());
  Hist1D::SingleHist1D* den_sh_nb0 = static_cast<Hist1D::SingleHist1D*>(den_h_nb0->GetComponent(proc.get()));
  TH1D den_nb0 = den_sh_nb0->scaled_hist_;
  Hist1D * num_h_nb0 = static_cast<Hist1D*>(pm->Figures()[first_index+7].get());
  Hist1D::SingleHist1D* num_sh_nb0 = static_cast<Hist1D::SingleHist1D*>(num_h_nb0->GetComponent(proc.get()));
  TH1D num_nb0 = num_sh_nb0->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_nb0 = new TGraphAsymmErrors(&num_nb0,&den_nb0,"cp");
  //nb1
  Hist1D * den_h_nb1 = static_cast<Hist1D*>(pm->Figures()[first_index+8].get());
  Hist1D::SingleHist1D* den_sh_nb1 = static_cast<Hist1D::SingleHist1D*>(den_h_nb1->GetComponent(proc.get()));
  TH1D den_nb1 = den_sh_nb1->scaled_hist_;
  Hist1D * num_h_nb1 = static_cast<Hist1D*>(pm->Figures()[first_index+9].get());
  Hist1D::SingleHist1D* num_sh_nb1 = static_cast<Hist1D::SingleHist1D*>(num_h_nb1->GetComponent(proc.get()));
  TH1D num_nb1 = num_sh_nb1->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_nb1 = new TGraphAsymmErrors(&num_nb1,&den_nb1,"cp");
  //nb2
  Hist1D * den_h_nb2 = static_cast<Hist1D*>(pm->Figures()[first_index+10].get());
  Hist1D::SingleHist1D* den_sh_nb2 = static_cast<Hist1D::SingleHist1D*>(den_h_nb2->GetComponent(proc.get()));
  TH1D den_nb2 = den_sh_nb2->scaled_hist_;
  Hist1D * num_h_nb2 = static_cast<Hist1D*>(pm->Figures()[first_index+11].get());
  Hist1D::SingleHist1D* num_sh_nb2 = static_cast<Hist1D::SingleHist1D*>(num_h_nb2->GetComponent(proc.get()));
  TH1D num_nb2 = num_sh_nb2->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_nb2 = new TGraphAsymmErrors(&num_nb2,&den_nb2,"cp");
  //drmax0
  Hist1D * den_h_drmax0 = static_cast<Hist1D*>(pm->Figures()[first_index+12].get());
  Hist1D::SingleHist1D* den_sh_drmax0 = static_cast<Hist1D::SingleHist1D*>(den_h_drmax0->GetComponent(proc.get()));
  TH1D den_drmax0 = den_sh_drmax0->scaled_hist_;
  Hist1D * num_h_drmax0 = static_cast<Hist1D*>(pm->Figures()[first_index+13].get());
  Hist1D::SingleHist1D* num_sh_drmax0 = static_cast<Hist1D::SingleHist1D*>(num_h_drmax0->GetComponent(proc.get()));
  TH1D num_drmax0 = num_sh_drmax0->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_drmax0 = new TGraphAsymmErrors(&num_drmax0,&den_drmax0,"cp");
  //drmax22
  Hist1D * den_h_drmax22 = static_cast<Hist1D*>(pm->Figures()[first_index+14].get());
  Hist1D::SingleHist1D* den_sh_drmax22 = static_cast<Hist1D::SingleHist1D*>(den_h_drmax22->GetComponent(proc.get()));
  TH1D den_drmax22 = den_sh_drmax22->scaled_hist_;
  Hist1D * num_h_drmax22 = static_cast<Hist1D*>(pm->Figures()[first_index+15].get());
  Hist1D::SingleHist1D* num_sh_drmax22 = static_cast<Hist1D::SingleHist1D*>(num_h_drmax22->GetComponent(proc.get()));
  TH1D num_drmax22 = num_sh_drmax22->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_drmax22 = new TGraphAsymmErrors(&num_drmax22,&den_drmax22,"cp");
  //fixht
  Hist1D * den_h_fixht = static_cast<Hist1D*>(pm->Figures()[first_index+16].get());
  Hist1D::SingleHist1D* den_sh_fixht = static_cast<Hist1D::SingleHist1D*>(den_h_fixht->GetComponent(proc.get()));
  TH1D den_fixht = den_sh_fixht->scaled_hist_;
  Hist1D * num_h_fixht = static_cast<Hist1D*>(pm->Figures()[first_index+17].get());
  Hist1D::SingleHist1D* num_sh_fixht = static_cast<Hist1D::SingleHist1D*>(num_h_fixht->GetComponent(proc.get()));
  TH1D num_fixht = num_sh_fixht->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_fixht = new TGraphAsymmErrors(&num_fixht,&den_fixht,"cp");
  //ambb0
  Hist1D * den_h_ambb0 = static_cast<Hist1D*>(pm->Figures()[first_index+18].get());
  Hist1D::SingleHist1D* den_sh_ambb0 = static_cast<Hist1D::SingleHist1D*>(den_h_ambb0->GetComponent(proc.get()));
  TH1D den_ambb0 = den_sh_ambb0->scaled_hist_;
  Hist1D * num_h_ambb0 = static_cast<Hist1D*>(pm->Figures()[first_index+19].get());
  Hist1D::SingleHist1D* num_sh_ambb0 = static_cast<Hist1D::SingleHist1D*>(num_h_ambb0->GetComponent(proc.get()));
  TH1D num_ambb0 = num_sh_ambb0->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_ambb0 = new TGraphAsymmErrors(&num_ambb0,&den_ambb0,"cp");
  //ambb140
  Hist1D * den_h_ambb140 = static_cast<Hist1D*>(pm->Figures()[first_index+20].get());
  Hist1D::SingleHist1D* den_sh_ambb140 = static_cast<Hist1D::SingleHist1D*>(den_h_ambb140->GetComponent(proc.get()));
  TH1D den_ambb140 = den_sh_ambb140->scaled_hist_;
  Hist1D * num_h_ambb140 = static_cast<Hist1D*>(pm->Figures()[first_index+21].get());
  Hist1D::SingleHist1D* num_sh_ambb140 = static_cast<Hist1D::SingleHist1D*>(num_h_ambb140->GetComponent(proc.get()));
  TH1D num_ambb140 = num_sh_ambb140->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_ambb140 = new TGraphAsymmErrors(&num_ambb140,&den_ambb140,"cp");
  //eltrig
  Hist1D * den_h_eltrig = static_cast<Hist1D*>(pm->Figures()[first_index+22].get());
  Hist1D::SingleHist1D* den_sh_eltrig = static_cast<Hist1D::SingleHist1D*>(den_h_eltrig->GetComponent(proc.get()));
  TH1D den_eltrig = den_sh_eltrig->scaled_hist_;
  Hist1D * num_h_eltrig = static_cast<Hist1D*>(pm->Figures()[first_index+23].get());
  Hist1D::SingleHist1D* num_sh_eltrig = static_cast<Hist1D::SingleHist1D*>(num_h_eltrig->GetComponent(proc.get()));
  TH1D num_eltrig = num_sh_eltrig->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_eltrig = new TGraphAsymmErrors(&num_eltrig,&den_eltrig,"cp");
  //jettrig
  Hist1D * den_h_jettrig = static_cast<Hist1D*>(pm->Figures()[first_index+24].get());
  Hist1D::SingleHist1D* den_sh_jettrig = static_cast<Hist1D::SingleHist1D*>(den_h_jettrig->GetComponent(proc.get()));
  TH1D den_jettrig = den_sh_jettrig->scaled_hist_;
  Hist1D * num_h_jettrig = static_cast<Hist1D*>(pm->Figures()[first_index+25].get());
  Hist1D::SingleHist1D* num_sh_jettrig = static_cast<Hist1D::SingleHist1D*>(num_h_jettrig->GetComponent(proc.get()));
  TH1D num_jettrig = num_sh_jettrig->scaled_hist_;
  TGraphAsymmErrors* plot_ratio_jettrig = new TGraphAsymmErrors(&num_jettrig,&den_jettrig,"cp");
  ////nel1
  //Hist1D * den_h_nel1 = static_cast<Hist1D*>(pm->Figures()[first_index+26].get());
  //Hist1D::SingleHist1D* den_sh_nel1 = static_cast<Hist1D::SingleHist1D*>(den_h_nel1->GetComponent(proc.get()));
  //TH1D den_nel1 = den_sh_nel1->scaled_hist_;
  //Hist1D * num_h_nel1 = static_cast<Hist1D*>(pm->Figures()[first_index+27].get());
  //Hist1D::SingleHist1D* num_sh_nel1 = static_cast<Hist1D::SingleHist1D*>(num_h_nel1->GetComponent(proc.get()));
  //TH1D num_nel1 = num_sh_nel1->scaled_hist_;
  //TGraphAsymmErrors* plot_ratio_nel1 = new TGraphAsymmErrors(&num_nel1,&den_nel1,"cp");
  ////nel0
  //Hist1D * den_h_nel0 = static_cast<Hist1D*>(pm->Figures()[first_index+28].get());
  //Hist1D::SingleHist1D* den_sh_nel0 = static_cast<Hist1D::SingleHist1D*>(den_h_nel0->GetComponent(proc.get()));
  //TH1D den_nel0 = den_sh_nel0->scaled_hist_;
  //Hist1D * num_h_nel0 = static_cast<Hist1D*>(pm->Figures()[first_index+29].get());
  //Hist1D::SingleHist1D* num_sh_nel0 = static_cast<Hist1D::SingleHist1D*>(num_h_nel0->GetComponent(proc.get()));
  //TH1D num_nel0 = num_sh_nel0->scaled_hist_;
  //TGraphAsymmErrors* plot_ratio_nel0 = new TGraphAsymmErrors(&num_nel0,&den_nel0,"cp");
  
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
  table_file.open("tables/"+out_filename+"_trig_sys_"+std::to_string(year)+".tex");
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
  table_file << make_systematic_table_row(plot_ratio_drmax0,"$\\Delta R_\\text{max}< 2.2$") << "\n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_drmax0));
  table_file << make_systematic_table_row(plot_ratio_drmax22,"$\\Delta R_\\text{max}\\geq 2.2$") << " \\hline \n";
  variations_table.push_back(get_systematic_table_values(plot_ratio_drmax22));
  variations_row = get_systematic_table_variations(variations_table);
  greater_variations_rows.push_back(variations_row);
  table_file << make_systematic_table_row(variations_row) << "\n";
  table_file << "\\multicolumn{8}{c}{$\\langle m_{bb}\\rangle$ variations (for $300 \\leq H_{T} < 400$~GeV)}\\\\ \\hline \n";
  table_file << make_systematic_table_row(plot_ratio_fixht,"$\\langle m_{bb}\\rangle$ Inclusive (Nominal)") << "\n";
  ambb_table.push_back(get_systematic_table_values(plot_ratio_fixht));
  table_file << make_systematic_table_row(plot_ratio_ambb0,"$\\langle m_{bb}\\rangle \\leq$ 140 GeV") << "\n";
  ambb_table.push_back(get_systematic_table_values(plot_ratio_ambb0));
  table_file << make_systematic_table_row(plot_ratio_ambb140,"$\\langle m_{bb}\\rangle >$ 140 GeV") << " \\hline \n";
  ambb_table.push_back(get_systematic_table_values(plot_ratio_ambb140));
  ambb_row = get_systematic_table_variations(ambb_table);
  greater_variations_rows.push_back(ambb_row);
  table_file << make_systematic_table_row(ambb_row) << "\n";
  //table_file << "\\multicolumn{8}{c}{$N_\\text{l}$ variations (Jet500 trigger, Max $p_\\text{Tjet}>500$~GeV)}\\\\ \\hline \n";
  //table_file << make_systematic_table_row(plot_ratio_nel1,"$N_\\text{e}=1$ (Nominal)") << "\n";
  //nel_table.push_back(get_systematic_table_values(plot_ratio_nel1));
  //table_file << make_systematic_table_row(plot_ratio_nel0,"$N_\\text{\\mu}=1$") << " \\hline \n";
  //nel_table.push_back(get_systematic_table_values(plot_ratio_nel0));
  //nel_row = get_systematic_table_variations(nel_table);
  //greater_variations_rows.push_back(nel_row);
  //table_file << make_systematic_table_row(nel_row) << "\n";
  systematics_rows.push_back(max_variations(greater_variations_rows));
  table_file << make_systematic_table_row(max_variations(greater_variations_rows)) << "\n";
  table_file << "\\multicolumn{8}{c}{Reference trigger (for Max $p_\\text{Tjet}>500$~GeV)}\\\\ \\hline \n";
  table_file << make_systematic_table_row(plot_ratio_eltrig,"SingleElectron (Nominal)") << "\n";
  reference_table.push_back(get_systematic_table_values(plot_ratio_eltrig));
  table_file << make_systematic_table_row(plot_ratio_jettrig,"Jet500") << " \\hline \n";
  reference_table.push_back(get_systematic_table_values(plot_ratio_jettrig));
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
  std::cout << "Wrote tables/higgsino/"+out_filename+"_"+std::to_string(year)+"_trig_sys.tex" << std::endl;
  return ret;
}

std::vector<std::vector<double>> generate_syst_extrapolation_table(int year, PlotMaker* pm, int first_index, shared_ptr<Process> proc) {
  //generate extrapolation systematics
  int pm_idx = first_index;
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
      EfficiencyPlot * eff_plot = static_cast<EfficiencyPlot*>(pm->Figures()[pm_idx].get());
      EfficiencyPlot::SingleEfficiencyPlot* singleff_plot = static_cast<EfficiencyPlot::SingleEfficiencyPlot*>(eff_plot->GetComponent(proc.get()));
      TH1D den_hist = singleff_plot->raw_denominator_hist_;
      TH1D num_hist = singleff_plot->raw_numerator_hist_;
      std::cout << "DEBUG: fitting " << pm_idx << std::endl;
      TFitResultPtr frp = eff_plot->data_ratio_plots_[0]->Fit(&linear,"S","",20,50);
      //std::cout << "Fit result, ht bin " << ht_bin_idx << ", met bin " << met_bin_idx << ": " << frp->Parameter(0) << ", " << frp->Parameter(1) << std::endl;
      eff_plot->data_ratio_plots_[0]->Draw("AP");
      linear.Draw("same");
      c->Update();
      c->SaveAs(("plots/trig_eff_0l_fit_htbin"+std::to_string(ht_bin_idx)+"_metbin"+std::to_string(met_bin_idx)+"_"+year_string+".pdf").c_str());
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
      pm_idx += 1;
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
  table_file.open("tables/"+out_filename+"_trig_syst_extrapolation_"+std::to_string(year)+".tex");
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
  std::cout << "Wrote tables/"+out_filename+"_trig_syst_extrapolation_"+std::to_string(year)+".tex" << std::endl;

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
void trig_postprocess_efficiencies(PlotMaker* pm, int first_index, shared_ptr<Process> met_proc, shared_ptr<Process> lep_proc, std::string year, std::string img_file_extension, std::vector<double> systematics, std::vector<std::vector<double>> extrapolation_systematics) {
  std::string lumi_string = "35.9 fb^{-1} (13 TeV)";
  if (year == "2017") {
  	lumi_string = "41.5 fb^{-1} (13 TeV)";
  }
  else if (year == "2018") {
  	lumi_string = "60 fb^{-1} (13 TeV)";
  }
  std::string met_trigger_string = "MET Triggers";
  std::string el_trigger_string = "Single Electron Triggers";
  std::string mu_trigger_string = "Single Muon Triggers";
  ofstream apply_effs_file;
  apply_effs_file.open(("src/higgsino/apply_trigeffs_"+out_filename+".cpp.txt").c_str());
  apply_effs_file << "#include <vector>\n";
  apply_effs_file << "#include \"core/baby.hpp\"\n";
  apply_effs_file << "#include \"core/process.hpp\"\n";
  apply_effs_file << "#include \"core/named_func.hpp\"\n";
  apply_effs_file << "#include \"higgsino/hig_utilities.hpp\"\n";
  //apply_effs_file << "#include \"higgsino/hig_functions.hpp\"\n";
  apply_effs_file << "#include \"higgsino/apply_trigeffs" << year << ".hpp\"\n\n";
  apply_effs_file << "namespace Higfuncs{\n\n";
  //make graphs
  //true met histograms
  //draw_2d_efficiency_graphs(pm, first_index, first_index+1, ht_bins, true_met_bins, "ht", "H_{T}", "b.ht()", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{SingleElectron}, N_{j}#geq 3, high #Delta #phi, N_{e}=1", lumi_string, met_trigger_string, "truemet_effhtbin", apply_effs_file, "get_0l_trigeff", year, img_file_extension, systematics);
  draw_2d_efficiency_graphs_variable_binning(pm, {first_index, first_index+2, first_index+4, first_index+6, first_index+8, first_index+10}, {first_index+1, first_index+3, first_index+5, first_index+7, first_index+9, first_index+11}, met_proc, ht_true_met_bins, "ht", "H_{T}", "b.ht()", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{SingleElectron}, N_{j}#geq 3, high #Delta #phi, N_{e}=1, p_{Te}<35 GeV", lumi_string, met_trigger_string, "truemet_effhtbin", apply_effs_file, "get_0l_trigeff", year, img_file_extension, systematics, extrapolation_systematics);
  //fake met histograms
  draw_2d_efficiency_graphs_variable_binning(pm, {first_index+12,first_index+14,first_index+16,first_index+18,first_index+20,first_index+22,first_index+24}, {first_index+13,first_index+15,first_index+17,first_index+19,first_index+21,first_index+23,first_index+25}, met_proc, ht_fake_met_bins, "ht", "H_{T}", "b.ht()", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{JetHT}, low #Delta #phi, N_{l veto}=0", lumi_string, met_trigger_string, "fakemet_effhtbin", apply_effs_file, "get_0l_fakemet_trigeff", year, img_file_extension, {}, {});
  if (do_controlregions) {
    //1/2l CR graphs
    draw_3d_efficiency_graphs(pm, {first_index+26, first_index+28, first_index+30}, {first_index+27, first_index+29, first_index+31}, singlelep_ht_bins, el_pt_bins, twodim_met_bins, "ht", "H_{T}", "b.ht()", "el_pt", "p_{Te}", "HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig())", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{JetHT}, N_{j}#geq 2, N_{e}=1", lumi_string, met_trigger_string+"|"+el_trigger_string, "elmet_effelptbin", apply_effs_file, "get_1el_trigeff", year, img_file_extension);
    draw_3d_efficiency_graphs(pm, {first_index+32, first_index+34, first_index+36}, {first_index+33, first_index+35, first_index+37}, singlelep_ht_bins, mu_pt_bins, twodim_met_bins, "ht", "H_{T}", "b.ht()", "mu_pt", "p_{T#mu}", "HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig())", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{JetHT}, N_{j}#geq 2, N_{#mu}=1", lumi_string, met_trigger_string+"|"+mu_trigger_string, "mumet_effmuptbin", apply_effs_file, "get_1mu_trigeff", year, img_file_extension);
    draw_1d_efficiency_graphs(pm, first_index+38, first_index+39, lep_proc, twoel_pt_bins, "el_pt", "Max p_{Te}", "HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig())","Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{e}=2",lumi_string,el_trigger_string,"elel_eff",apply_effs_file,"get_2el_trigeff",year,img_file_extension);
    draw_1d_efficiency_graphs(pm, first_index+40, first_index+41, lep_proc, twomu_pt_bins, "mu_pt", "Max p_{T#mu}", "HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig())","Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{#mu}=2",lumi_string,mu_trigger_string,"mumu_eff",apply_effs_file,"get_2mu_trigeff",year,img_file_extension);
  }
  apply_effs_file << "}\n";
  apply_effs_file.close();
  std::cout << "Wrote src/higgsino/apply_trigeffs_"+out_filename+".cpp" << std::endl;
}


void draw_3d_efficiency_graphs(PlotMaker* pm, std::vector<int> den_idx, std::vector<int> num_idx, std::vector<double> z_bins, std::vector<double> y_bins, std::vector<double> x_bins, std::string z_var_name, std::string z_title_name, std::string z_expression, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension) {
  //loop over ht - each ht bin is a new TGraphAsymmErrors
  std::vector<std::vector<TGraphAsymmErrors*>> new_eff_plots;
  for (unsigned int z_bin = 0; z_bin < z_bins.size()-1; z_bin++) {
    std::vector<TGraphAsymmErrors*> new_2d_eff_plots;
    Hist2D * effhist_den = static_cast<Hist2D*>(pm->Figures()[den_idx[z_bin]].get());
    TH2D denominator_hist = effhist_den->GetDataHist();
    Hist2D * effhist_num = static_cast<Hist2D*>(pm->Figures()[num_idx[z_bin]].get());
    TH2D numerator_hist = effhist_num->GetDataHist();
    for (unsigned int y_bin = 0; y_bin < y_bins.size()-1; y_bin++) {
      TH1D* numerator_1d_hist = numerator_hist.ProjectionX("numerator_1d_hist",y_bin+1,y_bin+1);
      TH1D* denominator_1d_hist = denominator_hist.ProjectionX("denominator_1d_hist",y_bin+1,y_bin+1);
      new_2d_eff_plots.push_back(new TGraphAsymmErrors(numerator_1d_hist,denominator_1d_hist));
      make_efficiency_plot(numerator_1d_hist,denominator_1d_hist,new_2d_eff_plots[y_bin],lumi_string,cuts_string+", "+y_title_name+"#in ["+std::to_string(y_bins[y_bin])+", "+std::to_string(y_bins[y_bin+1])+"], "+z_title_name+"#in ["+std::to_string(z_bins[z_bin])+", "+std::to_string(z_bins[z_bin+1])+"]",x_title_name,"Efficiency "+new_trigger_string,"GeV",plot_name+"_"+std::to_string(z_bin)+"_"+std::to_string(y_bin)+"_"+year+"."+img_file_extension,true,false);
    }
    new_eff_plots.push_back(new_2d_eff_plots);
  }
  //write out function to file
  out_file << "const NamedFunc " << func_name << year << "(\"" << func_name << year << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.; // Not used, but for reference\n";
  out_file << "  float eff = 1., " << x_var_name << " = " << x_expression << ", " << y_var_name << " = " << y_expression << ", " << z_var_name << " = " << z_expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  for (unsigned int x_bin = 0; x_bin < x_bins.size()-1; x_bin++) {
    double x_upper = x_bin==(x_bins.size()-2) ? 9999 : x_bins[x_bin+1];
    for (unsigned int y_bin = 0; y_bin < y_bins.size()-1; y_bin++) {
      double y_upper = y_bin==(y_bins.size()-2) ? 9999 : y_bins[y_bin+1];
      for (unsigned int z_bin = 0; z_bin < z_bins.size()-1; z_bin++) {
        double z_upper = z_bin==(z_bins.size()-2) ? 9999 : z_bins[z_bin+1];
        double* eff_pts = new_eff_plots[z_bin][y_bin]->GetY();
        if (first) {
          out_file << "  if (" << z_var_name << "> ";
          first = false;
        }
        else out_file << "  else if (" << z_var_name << "> ";
        out_file << z_bins[z_bin] << " && " << z_var_name << "<= " << z_upper;
        out_file << " && " << y_var_name << "> " << y_bins[y_bin] << " && " << y_var_name << "<= " << y_upper;
        out_file << " && " << x_var_name << "> " << x_bins[x_bin] << " && " << x_var_name << "<= " << x_upper;
        float errup = new_eff_plots[z_bin][y_bin]->GetErrorYhigh(x_bin);
        float errdown = new_eff_plots[z_bin][y_bin]->GetErrorYlow(x_bin);
        out_file << ") {eff = " << eff_pts[x_bin] << "; errup = " << errup << "; errdown = " << errdown << ";}\n";
      }
    }
  }
  out_file << "  std::vector<double> ret = {eff, errup, errdown};\n";
  out_file << "  return ret;\n";
  out_file << "});\n\n";
}


void draw_2d_efficiency_graphs(PlotMaker* pm, int den_idx, int num_idx, std::vector<double> y_bins, std::vector<double> x_bins, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension, std::vector<double> systematics) {
  //loop over ht - each ht bin is a new TGraphAsymmErrors
  std::vector<TGraphAsymmErrors*> new_eff_plots;
  Hist2D * effhist_den = static_cast<Hist2D*>(pm->Figures()[den_idx].get());
  TH2D denominator_hist = effhist_den->GetDataHist();
  Hist2D * effhist_num = static_cast<Hist2D*>(pm->Figures()[num_idx].get());
  TH2D numerator_hist = effhist_num->GetDataHist();
  for (unsigned int y_bin = 0; y_bin < y_bins.size()-1; y_bin++) {
    TH1D* numerator_1d_hist = numerator_hist.ProjectionX("numerator_1d_hist",y_bin+1,y_bin+1);
    TH1D* denominator_1d_hist = denominator_hist.ProjectionX("denominator_1d_hist",y_bin+1,y_bin+1);
    new_eff_plots.push_back(new TGraphAsymmErrors(numerator_1d_hist,denominator_1d_hist));
    make_efficiency_plot(numerator_1d_hist,denominator_1d_hist,new_eff_plots[y_bin],lumi_string,cuts_string+", "+y_title_name+"#in ["+std::to_string(y_bins[y_bin])+", "+std::to_string(y_bins[y_bin+1])+"]",x_title_name,"Efficiency "+new_trigger_string,"GeV",plot_name+"_"+std::to_string(y_bin)+"_"+year+"."+img_file_extension,true,false);
  }
  //write out function to file
  out_file << "const NamedFunc " << func_name << year << "(\"" << func_name << year << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.; // Not used, but for reference\n";
  out_file << "  float eff = 1., " << x_var_name << " = " << x_expression << ", " << y_var_name << " = " << y_expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  for (unsigned int x_bin = 0; x_bin < x_bins.size()-1; x_bin++) {
    double x_upper = x_bin==(x_bins.size()-2) ? 9999 : x_bins[x_bin+1];
    for (unsigned int y_bin = 0; y_bin < y_bins.size()-1; y_bin++) {
      double y_upper = y_bin==(y_bins.size()-2) ? 9999 : y_bins[y_bin+1];
      double* eff_pts = new_eff_plots[y_bin]->GetY();
      if (first) {
        out_file << "  if (" << y_var_name << "> ";
        first = false;
      }
      else out_file << "  else if (" << y_var_name << "> ";
      out_file << y_bins[y_bin] << " && " << y_var_name << "<= " << y_upper;
      out_file << " && " << x_var_name << "> " << x_bins[x_bin] << " && " << x_var_name << "<= " << x_upper;
      float errup = new_eff_plots[y_bin]->GetErrorYhigh(x_bin);
      float errdown = new_eff_plots[y_bin]->GetErrorYlow(x_bin);
      //figure out systematics bin and add errors in quadrature, making sure not to exceed 1 or 0
      if (systematics.size()>0) {
        for (unsigned int sys_bin = 0; sys_bin < (sys_met_bins.size()-1); sys_bin++) {
          if (sys_met_bins[sys_bin] <= x_bins[x_bin] && x_bins[x_bin] < sys_met_bins[sys_bin+1]) {
      	    errup = TMath::Sqrt(errup*errup+eff_pts[x_bin]*eff_pts[x_bin]*systematics[sys_bin]*systematics[sys_bin]);
      	    errdown = TMath::Sqrt(errdown*errdown+eff_pts[x_bin]*eff_pts[x_bin]*systematics[sys_bin]*systematics[sys_bin]);
      	    if (errdown < new_eff_plots[y_bin]->GetErrorYlow(x_bin)) {
      	      std::cout << "stats: " << new_eff_plots[y_bin]->GetErrorYlow(x_bin) << ", value: " << eff_pts[x_bin] << ", sys: " << systematics[sys_bin] << std::endl;
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
}

void draw_2d_efficiency_graphs_variable_binning(PlotMaker* pm, std::vector<int> den_idx, std::vector<int> num_idx, shared_ptr<Process> proc, std::vector<double> y_bins, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension, std::vector<double> systematics, std::vector<std::vector<double>> extrapolation_systematics) {
  //loop over ht - each ht bin is a new TGraphAsymmErrors
  std::vector<TGraphAsymmErrors*> new_eff_plots;
  std::vector<std::vector<double>> x_bins;
  for (unsigned int y_bin = 0; y_bin < y_bins.size()-1; y_bin++) {
    Hist1D * den_h1d = static_cast<Hist1D*>(pm->Figures()[den_idx[y_bin]].get());
    Hist1D::SingleHist1D* den_sh = static_cast<Hist1D::SingleHist1D*>(den_h1d->GetComponent(proc.get()));
    TH1D denominator_1d_hist = den_sh->scaled_hist_;
    Hist1D * num_h1d = static_cast<Hist1D*>(pm->Figures()[num_idx[y_bin]].get());
    Hist1D::SingleHist1D* num_sh = static_cast<Hist1D::SingleHist1D*>(num_h1d->GetComponent(proc.get()));
    TH1D numerator_1d_hist = num_sh->scaled_hist_;
    new_eff_plots.push_back(new TGraphAsymmErrors(&numerator_1d_hist,&denominator_1d_hist,"cp"));
    std::vector<double> x_bins_1d;
    for (unsigned int x_bin = 1; x_bin <= static_cast<unsigned int>(numerator_1d_hist.GetNbinsX()); x_bin++) {
      x_bins_1d.push_back(numerator_1d_hist.GetBinLowEdge(x_bin));
    }
    x_bins_1d.push_back(numerator_1d_hist.GetBinLowEdge(numerator_1d_hist.GetNbinsX())+numerator_1d_hist.GetBinWidth(numerator_1d_hist.GetNbinsX()));
    x_bins.push_back(x_bins_1d);
    make_efficiency_plot(&numerator_1d_hist,&denominator_1d_hist,new_eff_plots[y_bin],lumi_string,cuts_string+", "+y_title_name+"#in ["+std::to_string(y_bins[y_bin])+", "+std::to_string(y_bins[y_bin+1])+"]",x_title_name,"Efficiency "+new_trigger_string,"GeV",plot_name+"_"+std::to_string(y_bin)+"_"+year+"."+img_file_extension,true,false);
  }
  //write out function to file
  out_file << "const NamedFunc " << func_name << year << "(\"" << func_name << year << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.; // Not used, but for reference\n";
  out_file << "  float eff = 1., " << x_var_name << " = " << x_expression << ", " << y_var_name << " = " << y_expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  for (unsigned int y_bin = 0; y_bin < y_bins.size()-1; y_bin++) {
    double y_upper = y_bin==(y_bins.size()-2) ? 9999 : y_bins[y_bin+1];
    double* eff_pts = new_eff_plots[y_bin]->GetY();
    for (unsigned int x_bin = 0; x_bin < x_bins[y_bin].size()-1; x_bin++) {
      double x_upper = x_bin==(x_bins[y_bin].size()-2) ? 9999 : x_bins[y_bin][x_bin+1];
      if (first) {
        out_file << "  if (" << y_var_name << "> ";
        first = false;
      }
      else out_file << "  else if (" << y_var_name << "> ";
      out_file << y_bins[y_bin] << " && " << y_var_name << "<= " << y_upper;
      out_file << " && " << x_var_name << "> " << x_bins[y_bin][x_bin] << " && " << x_var_name << "<= " << x_upper;
      float errup = new_eff_plots[y_bin]->GetErrorYhigh(x_bin);
      float errdown = new_eff_plots[y_bin]->GetErrorYlow(x_bin);
      //figure out systematics bin and add errors in quadrature, making sure not to exceed 1 or 0
      if (systematics.size()>0) {
        for (unsigned int sys_bin = 0; sys_bin < (sys_met_bins.size()-1); sys_bin++) {
          if (sys_met_bins[sys_bin] <= x_bins[y_bin][x_bin] && x_bins[y_bin][x_bin] < sys_met_bins[sys_bin+1]) {
      	    errup = TMath::Sqrt(errup*errup+eff_pts[x_bin]*eff_pts[x_bin]*systematics[sys_bin]*systematics[sys_bin]);
      	    errdown = TMath::Sqrt(errdown*errdown+eff_pts[x_bin]*eff_pts[x_bin]*systematics[sys_bin]*systematics[sys_bin]);
      	    if (errdown < new_eff_plots[y_bin]->GetErrorYlow(x_bin)) {
      	      std::cout << "stats: " << new_eff_plots[y_bin]->GetErrorYlow(x_bin) << ", value: " << eff_pts[x_bin] << ", sys: " << systematics[sys_bin] << std::endl;
      	    }
      	  }
        }
        //apply extrapolation systematic
        for (unsigned int sys_ht_bin = 0; sys_ht_bin < (sys_ht_bins.size()-1); sys_ht_bin++) {
          for (unsigned int sys_met_bin = 0; sys_met_bin < (sys_met_bins.size()-1); sys_met_bin++) {
            if (sys_ht_bins[sys_ht_bin] <= y_bins[y_bin] && y_bins[y_bin] < sys_ht_bins[sys_ht_bin+1]) {
              if (sys_met_bins[sys_met_bin] <= x_bins[y_bin][x_bin] && x_bins[y_bin][x_bin] < sys_met_bins[sys_met_bin+1]) {
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
}

void draw_1d_efficiency_graphs(PlotMaker* pm, int den_idx, int num_idx, shared_ptr<Process> proc, std::vector<double> x_bins, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension) {
  Hist1D * den_h1d = static_cast<Hist1D*>(pm->Figures()[den_idx].get());
  Hist1D::SingleHist1D* den_sh = static_cast<Hist1D::SingleHist1D*>(den_h1d->GetComponent(proc.get()));
  TH1D denominator_hist = den_sh->scaled_hist_;
  Hist1D * num_h1d = static_cast<Hist1D*>(pm->Figures()[num_idx].get());
  Hist1D::SingleHist1D* num_sh = static_cast<Hist1D::SingleHist1D*>(num_h1d->GetComponent(proc.get()));
  TH1D numerator_hist = num_sh->scaled_hist_;
  TGraphAsymmErrors* new_eff_hist = new TGraphAsymmErrors(&numerator_hist,&denominator_hist,"cp");
  make_efficiency_plot(&numerator_hist,&denominator_hist,new_eff_hist,lumi_string,cuts_string,x_title_name,"Efficiency "+new_trigger_string,"GeV",plot_name+"_"+year+"."+img_file_extension,true,false);
  //write out functions to file
  out_file << "const NamedFunc " << func_name << year << "(\"" << func_name << year << "\", [](const Baby &b) -> NamedFunc::VectorType{\n";
  out_file << "  float errup=0., errdown=0.; // Not used, but for reference\n";
  out_file << "  float eff = 1., " << x_var_name << " = " << x_expression << ";\n";
  out_file << "  errup+=errdown; //suppress unused warning\n";
  bool first = true;
  double* eff_x_vals = new_eff_hist->GetY();
  for (unsigned int x_bin = 0; x_bin < x_bins.size()-1; x_bin++) {
    double x_upper = x_bin==(x_bins.size()-2) ? 9999 : x_bins[x_bin+1];
    if (first) {
      out_file << "  if (" << x_var_name << "> ";
      first = false;
    }
    else out_file << "  else if (" << x_var_name << "> ";
    out_file << x_bins[x_bin] << " && " << x_var_name << "<= " << x_upper;
    out_file << ") {eff = " << eff_x_vals[x_bin] << "; errup = " << new_eff_hist->GetErrorYhigh(x_bin) << "; errdown = " << new_eff_hist->GetErrorYlow(x_bin) << ";}\n";
  }
  out_file << "  std::vector<double> ret = {eff, errup, errdown};\n";
  out_file << "  return ret;\n";
  out_file << "});\n\n";
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {"out", required_argument, 0, 0},
      {"year", required_argument, 0, 0},
      {"plots", no_argument, 0, 0},
      {"systematics", no_argument, 0, 0},
      {"efficiencies", no_argument, 0, 0},
      {"includecr", no_argument, 0, 0},
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
      } else if (optname == "plots") {
        do_variables = true;
      } else if (optname == "systematics") {
        do_systematics = true;
      } else if (optname == "efficiencies") {
        do_efficiency = true;
      } else if (optname == "includecr") {
        do_controlregions = true;
      } else if (optname == "out") {
        out_filename = optarg;
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
