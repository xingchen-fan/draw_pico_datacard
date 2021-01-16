#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TMath.h"
#include "TRandom.h"
#include "TError.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "higgsino/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  //argument options
  bool single_thread = false;
  bool do_variables = true;
  bool do_systematics = true;
  bool do_efficiency = true;
  std::string year_string = "2016";
  std::string out_filename = "triggereff";
  bool do_controlregions = true; //false will omit plots of 1l and 2l CRs, speeding up processing by about 10~20x

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
}

//plot helper function declarations
void generate_ratio_plots(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, const char* hist_title, string hist_name);
void generate_2d_efficiencies(PlotMaker* pm, int denominator_index, int numerator_index, const char* hist_title, string hist_name);
void make_ambb_plots(PlotMaker* pm, int den_index_0, int num_index_0, int den_index_100, int num_index_100, int den_index150, int num_index150, shared_ptr<Process> proc, std::string lumi_string, std::string description, std::string name);
void make_efficiency_plot_wrapper(PlotMaker* pm, int denominator_index, int numerator_index, shared_ptr<Process> proc, std::string lumi_string, std::string title, std::string xaxis, std::string yaxis, std::string units, std::string out_filename, bool var_width_bins=false, bool is_simulation=false);
void make_efficiency_plot(TH1D* hist_num, TH1D* hist_den, TGraphAsymmErrors* hist_ratio, std::string lumi_string, std::string title, std::string xaxis, std::string yaxis, std::string units, std::string out_filename, bool var_width_bins=false, bool is_simulation=false);

//systematic helper function declarations
std::vector<double> trig_postprocess_systematics(int year, PlotMaker* pm, int first_index, shared_ptr<Process> proc);
std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name);
std::vector<double> get_systematic_table_values(TGraphAsymmErrors *gr);
std::vector<double> get_systematic_table_variations(std::vector<std::vector<double>> sys_table);
std::string make_systematic_table_row(std::vector<double> values, bool total=false);
std::vector<double> add_variations_quadrature(std::vector<std::vector<double>> values);
std::vector<double> max_variations(std::vector<std::vector<double>> values);

//efficiency helper function declarations
void trig_postprocess_efficiencies(PlotMaker* pm, int first_index, shared_ptr<Process> met_proc, shared_ptr<Process> lep_proc, std::string year, std::string img_file_extension, std::vector<double> systematics);
void draw_3d_efficiency_graphs(PlotMaker* pm, std::vector<int> den_idx, std::vector<int> num_idx, std::vector<double> z_bins, std::vector<double> y_bins, std::vector<double> x_bins, std::string z_var_name, std::string z_title_name, std::string z_expression, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension);
void draw_2d_efficiency_graphs(PlotMaker* pm, int den_idx, int num_idx, std::vector<double> y_bins, std::vector<double> x_bins, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension, std::vector<double> systematics);
void draw_2d_efficiency_graphs_variable_binning(PlotMaker* pm, std::vector<int> den_idx, std::vector<int> num_idx, shared_ptr<Process> proc, std::vector<double> y_bins, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension);
void draw_1d_efficiency_graphs(PlotMaker* pm, int den_idx, int num_idx, shared_ptr<Process> proc, std::vector<double> x_bins, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  // Define 1D+2D plot types of interest
  Palette colors("txt/colors.txt", "default");
  PlotOpt lin_lumi("txt/plot_styles.txt", "Std1D");
  lin_lumi.Title(TitleType::info).Overflow(OverflowType::overflow);
  vector<PlotOpt> all_plot_types = {lin_lumi};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style2D().Title(TitleType::info).Overflow(OverflowType::overflow)};

  //year-based settings, default to 2016 but changed below
  double lumi = 35.9;
  int year = 2016;
  string data_dir = "/net/cms25/cms25r0/pico/NanoAODv7/higgsino_inyo/2016/data/";
  string mc_dir = "/net/cms29/cms29r0/pico/NanoAODv5/higgsino_humboldt/2016/mc/";
  string singleelectron_name = "SingleElectron";

  if (year_string=="2016") year = 2016;
  else if (year_string=="2017") year = 2017;
  else if (year_string=="2018") year = 2018;
  else {
    std::cout << "ERROR: unsupported year." << std::endl;
    return 1;
  }

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
  //set<string> str_met150({data_dir+"skim_met150/raw_pico_met150_MET__Run2016C*.root"});
  pro_met150 = Process::MakeShared<Baby_pico>("Data MET150", Process::Type::data, kBlack, str_met150, "stitch");
  set<string> str_met150_mc({mc_dir+"skim_met150/pico_met150_TTJets_SingleLeptFromT*.root"});
  pro_met150_mc = Process::MakeShared<Baby_pico>("TTbar 1l", Process::Type::background, kBlack, str_met150_mc, "stitch");
  set<string> str_1l2j({data_dir+"skim_1l2j/raw_pico_1l2j_"+singleelectron_name+"*.root",data_dir+"skim_1l2j/raw_pico_1l2j_MET*.root",data_dir+"skim_1l2j/raw_pico_1l2j_SingleMuon*.root",data_dir+"skim_1l2j/raw_pico_1l2j_JetHT*.root"});
  pro_1l2j = Process::MakeShared<Baby_pico>("2016 Data", Process::Type::data, kBlack, str_1l2j, "stitch");
  procs_met150.push_back(pro_met150);
  procs_met150_mc.push_back(pro_met150_mc);
  procs_1l2j.push_back(pro_1l2j);

  //named funcs
  //const NamedFunc met_trigger = "HLT_PFMET110_PFMHT110_IDTight||HLT_PFMETNoMu110_PFMHTNoMu110_IDTight||HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight||HLT_PFMET120_PFMHT120_IDTight_PFHT60||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60";
  //const NamedFunc el_trigger = "HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf||HLT_Ele115_CaloIdVT_GsfTrkIdT";
  //const NamedFunc mu_trigger = "HLT_IsoMu24||HLT_IsoMu27||HLT_Mu50";
  //const NamedFunc ht_trigger = "HLT_PFJet500";

  //const NamedFunc met_trigger = "HLT_PFMET90_PFMHT90_IDTight||HLT_PFMETNoMu90_PFMHTNoMu90_IDTight||HLT_PFMET100_PFMHT100_IDTight||HLT_PFMETNoMu100_PFMHTNoMu100_IDTight||HLT_PFMET110_PFMHT110_IDTight||HLT_PFMETNoMu110_PFMHTNoMu110_IDTight||HLT_PFMET120_PFMHT120_IDTight||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight||HLT_PFMET130_PFMHT130_IDTight||HLT_PFMETNoMu130_PFMHTNoMu130_IDTight||HLT_PFMET140_PFMHT140_IDTight||HLT_PFMETNoMu140_PFMHTNoMu140_IDTight||HLT_PFMET100_PFMHT100_IDTight_PFHT60||HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60||HLT_PFMET110_PFMHT110_IDTight_PFHT60||HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60||HLT_PFMET120_PFMHT120_IDTight_PFHT60||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60||HLT_PFMET130_PFMHT130_IDTight_PFHT60||HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60||HLT_PFMET140_PFMHT140_IDTight_PFHT60||HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60||HLT_PFMET120_PFMHT120_IDTight_HFCleaned||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned||HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned";
  //const NamedFunc el_trigger = "HLT_Ele25_WPTight_Gsf||HLT_Ele27_WPTight_Gsf||HLT_Ele28_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele32_WPTight_Gsf||HLT_Ele35_WPTight_Gsf||HLT_Ele45_WPLoose_Gsf||HLT_Ele105_CaloIdVT_GsfTrkIdT||HLT_Ele115_CaloIdVT_GsfTrkIdT||HLT_Ele135_CaloIdVT_GsfTrkIdT||HLT_Ele145_CaloIdVT_GsfTrkIdT||HLT_Ele25_eta2p1_WPTight_Gsf||HLT_Ele27_eta2p1_WPTight_Gsf||HLT_Ele27_eta2p1_WPLoose_Gsf||HLT_Ele20_WPLoose_Gsf||HLT_Ele20_eta2p1_WPLoose_Gsf||HLT_Ele25_eta2p1_WPLoose_Gsf||HLT_Ele15_IsoVVVL_PFHT350||HLT_Ele15_IsoVVVL_PFHT400||HLT_Ele15_IsoVVVL_PFHT450||HLT_Ele15_IsoVVVL_PFHT600||HLT_Ele50_IsoVVVL_PFHT450";
  //const NamedFunc mu_trigger = "HLT_IsoMu20||HLT_IsoMu22||HLT_IsoMu24||HLT_IsoMu27||HLT_IsoTkMu20||HLT_IsoTkMu22||HLT_IsoTkMu24||HLT_Mu50||HLT_Mu55||HLT_TkMu50||HLT_IsoMu22_eta2p1||HLT_IsoMu24_eta2p1||HLT_Mu45_eta2p1||HLT_Mu15_IsoVVVL_PFHT350||HLT_Mu15_IsoVVVL_PFHT400||HLT_Mu15_IsoVVVL_PFHT450||HLT_Mu15_IsoVVVL_PFHT600||HLT_Mu50_IsoVVVL_PFHT400||HLT_Mu50_IsoVVVL_PFHT450";
  const NamedFunc jet_trigger = "HLT_PFJet500";
  //const NamedFunc ht_trigger = "(HLT_PFHT125 || HLT_PFHT200 || HLT_PFHT300 || HLT_PFHT400 || HLT_PFHT475 || HLT_PFHT600 || HLT_PFHT650 || HLT_PFHT800 || HLT_PFHT900 || HLT_PFHT180 || HLT_PFHT370 || HLT_PFHT430 || HLT_PFHT510 || HLT_PFHT590 || HLT_PFHT680 || HLT_PFHT780 || HLT_PFHT890 || HLT_PFHT1050 || HLT_PFHT250 || HLT_PFHT350)";

  const NamedFunc met_trigger("met_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_met_trigger = b.HLT_PFMET90_PFMHT90_IDTight()||b.HLT_PFMETNoMu90_PFMHTNoMu90_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight()||b.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight()||b.HLT_PFMET110_PFMHT110_IDTight()||b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight()||b.HLT_PFMET120_PFMHT120_IDTight()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight()||b.HLT_PFMET130_PFMHT130_IDTight()||b.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight()||b.HLT_PFMET140_PFMHT140_IDTight()||b.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight_PFHT60()||b.HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60()||b.HLT_PFMET110_PFMHT110_IDTight_PFHT60()||b.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60()||b.HLT_PFMET130_PFMHT130_IDTight_PFHT60()||b.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60()||b.HLT_PFMET140_PFMHT140_IDTight_PFHT60()||b.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_HFCleaned()||b.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned();
    return r_met_trigger;
  });

  const NamedFunc metnomu_trigger("metnomu_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_met_trigger = b.HLT_PFMET90_PFMHT90_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight()||b.HLT_PFMET110_PFMHT110_IDTight()||b.HLT_PFMET120_PFMHT120_IDTight()||b.HLT_PFMET130_PFMHT130_IDTight()||b.HLT_PFMET140_PFMHT140_IDTight()||b.HLT_PFMET100_PFMHT100_IDTight_PFHT60()||b.HLT_PFMET110_PFMHT110_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60()||b.HLT_PFMET130_PFMHT130_IDTight_PFHT60()||b.HLT_PFMET140_PFMHT140_IDTight_PFHT60()||b.HLT_PFMET120_PFMHT120_IDTight_HFCleaned()||b.HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned();
    return r_met_trigger;
  });

  const NamedFunc el_trigger("el_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_el_trigger = b.HLT_Ele25_WPTight_Gsf()||b.HLT_Ele27_WPTight_Gsf()||b.HLT_Ele28_WPTight_Gsf()||b.HLT_Ele32_WPTight_Gsf_L1DoubleEG()||b.HLT_Ele32_WPTight_Gsf()||b.HLT_Ele35_WPTight_Gsf()||b.HLT_Ele45_WPLoose_Gsf()||b.HLT_Ele105_CaloIdVT_GsfTrkIdT()||b.HLT_Ele115_CaloIdVT_GsfTrkIdT()||b.HLT_Ele135_CaloIdVT_GsfTrkIdT()||b.HLT_Ele145_CaloIdVT_GsfTrkIdT()||b.HLT_Ele25_eta2p1_WPTight_Gsf()||b.HLT_Ele27_eta2p1_WPTight_Gsf()||b.HLT_Ele27_eta2p1_WPLoose_Gsf()||b.HLT_Ele20_WPLoose_Gsf()||b.HLT_Ele20_eta2p1_WPLoose_Gsf()||b.HLT_Ele25_eta2p1_WPLoose_Gsf()||b.HLT_Ele15_IsoVVVL_PFHT350()||b.HLT_Ele15_IsoVVVL_PFHT400()||b.HLT_Ele15_IsoVVVL_PFHT450()||b.HLT_Ele15_IsoVVVL_PFHT600()||b.HLT_Ele50_IsoVVVL_PFHT450();
    return r_el_trigger;
  });

  const NamedFunc mu_trigger("mu_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_mu_trigger = b.HLT_IsoMu20()||b.HLT_IsoMu22()||b.HLT_IsoMu24()||b.HLT_IsoMu27()||b.HLT_IsoTkMu20()||b.HLT_IsoTkMu22()||b.HLT_IsoTkMu24()||b.HLT_Mu50()||b.HLT_Mu55()||b.HLT_TkMu50()||b.HLT_IsoMu22_eta2p1()||b.HLT_IsoMu24_eta2p1()||b.HLT_Mu45_eta2p1()||b.HLT_Mu15_IsoVVVL_PFHT350()||b.HLT_Mu15_IsoVVVL_PFHT400()||b.HLT_Mu15_IsoVVVL_PFHT450()||b.HLT_Mu15_IsoVVVL_PFHT600()||b.HLT_Mu50_IsoVVVL_PFHT400()||b.HLT_Mu50_IsoVVVL_PFHT450();
    return r_mu_trigger;
  });

  const NamedFunc ht_trigger("ht_trigger", [](const Baby &b) -> NamedFunc::ScalarType{
    //can't used string-based named func because name being too long causes histograms to crash
    bool r_ht_trigger = b.HLT_PFHT125()||b.HLT_PFHT200()||b.HLT_PFHT300()||b.HLT_PFHT400()||b.HLT_PFHT475()||b.HLT_PFHT600()||b.HLT_PFHT650()||b.HLT_PFHT800()||b.HLT_PFHT900()||b.HLT_PFHT180()||b.HLT_PFHT370()||b.HLT_PFHT430()||b.HLT_PFHT510()||b.HLT_PFHT590()||b.HLT_PFHT680()||b.HLT_PFHT780()||b.HLT_PFHT890()||b.HLT_PFHT1050()||b.HLT_PFHT250()||b.HLT_PFHT350();
    return r_ht_trigger;
  });

  const NamedFunc jetht_trigger = "HLT_PFJet500" || ht_trigger;

  const NamedFunc ht_loose_jets("ht_loose_jets", [](const Baby &b) -> NamedFunc::ScalarType{
		  float r_ht_loose_jets = 0;
		  for (unsigned jet_idx = 0; jet_idx < b.jet_pt()->size(); jet_idx++) {
		    //jetid is bugged in current nano2pico
		    if (b.jet_isgood()->at(jet_idx) && (b.jet_id()->at(jet_idx))) {
		      r_ht_loose_jets += b.jet_pt()->at(jet_idx);
		    }
		  }
		  return r_ht_loose_jets;
  });

  const NamedFunc high_pt_jet("high_pt_jet", [](const Baby &b) -> NamedFunc::ScalarType{
		  bool r_high_pt_jet = false;
		  if (b.njet()>0) {
		    if (b.jet_pt()->at(0) > 500) {
			    r_high_pt_jet = true;
		    }
		  }
		  return r_high_pt_jet;
  });

  const NamedFunc nbb("nbb", [](const Baby &b) -> NamedFunc::ScalarType{
      int r_nbb = 0;
      for (unsigned int fjet_idx = 0; fjet_idx < static_cast<unsigned int>(b.nfjet()); fjet_idx++) {
        if (fjet_idx >= 2)  break;
        if (b.fjet_deep_md_hbb_btv()->at(fjet_idx)>0.7)
          r_nbb++;
      }
		  return r_nbb;
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

  const NamedFunc reminiaod_ht("reminiaod_ht", [](const Baby &b) -> NamedFunc::ScalarType{
    float pico_ht = b.ht();
    TRandom3 rndgen;
    float deltaht_mean = 7.67, deltaht_sigma = 0;
    //float deltaht_mean = 7.67, deltaht_sigma = 18.79;
    //if (pico_ht> 0 && pico_ht <= 200) { deltaht_mean = -0.357143; deltaht_sigma = 18.681;}
    //else if (pico_ht> 200 && pico_ht <= 300) { deltaht_mean = 5.41667; deltaht_sigma = 9.86013;}
    //else if (pico_ht> 300 && pico_ht <= 400) { deltaht_mean = 5.18519; deltaht_sigma = 13.5985;}
    //else if (pico_ht> 400 && pico_ht <= 500) { deltaht_mean = 9.25; deltaht_sigma = 15.5141;}
    //else if (pico_ht> 500 && pico_ht <= 600) { deltaht_mean = 12.2143; deltaht_sigma = 13.0353;}
    //else if (pico_ht> 600 && pico_ht <= 9999) { deltaht_mean = 9.85294; deltaht_sigma = 13.7859;}
    return pico_ht-rndgen.Gaus(deltaht_mean, deltaht_sigma);
  });
  
  const NamedFunc reminiaod_met("reminiaod_met", [](const Baby &b) -> NamedFunc::ScalarType{
    float pico_met = b.met();
    TRandom3 rndgen;
    float deltamet_mean = 3., deltamet_sigma = 0.;
    //float deltamet_mean = 3.168, deltamet_sigma = 12.51;
    //if (pico_met> 0 && pico_met <= 80) { deltamet_mean = -1.36364; deltamet_sigma = 12.0226;}
    //else if (pico_met> 80 && pico_met <= 160) { deltamet_mean = 2.23333; deltamet_sigma = 9.4361;}
    //else if (pico_met> 160 && pico_met <= 240) { deltamet_mean = 4.66495; deltamet_sigma = 12.3294;}
    //else if (pico_met> 240 && pico_met <= 320) { deltamet_mean = 3.58108; deltamet_sigma = 13.0558;}
    //else if (pico_met> 320 && pico_met <= 9999) { deltamet_mean = 3.33333; deltamet_sigma = 1.86339;}
    return pico_met-rndgen.Gaus(deltamet_mean, deltamet_sigma);
  });

  const NamedFunc pass_filters("pass_filters", [](const Baby &b) -> NamedFunc::ScalarType{
    if (!b.pass_goodv() || !b.pass_hbhe() || !b.pass_hbheiso() || !b.pass_ecaldeadcell() || !b.pass_badpfmu() || !b.pass_muon_jet()) return false;
    if (b.type()/1000 == 0 && !b.pass_eebadsc()) return false; //only apply eebadsc fiter for data
    if ((b.type()/1000 != 106)  && !b.pass_cschalo_tight()) return false; //not for fastsim
    if (!b.pass_low_neutral_jet()) return false;
    if (!b.pass_htratio_dphi_tight()) return false;
    if (!b.pass_jets()) return false; //was modified
    if ((abs(b.SampleType())==2017 || abs(b.SampleType())==2018) && !Higfuncs::pass_ecalnoisejet.GetScalar(b)) return false; 
    if (!Higfuncs::pass_hemveto.GetScalar(b)) return false;
    return true;
  });

  PlotMaker pm;

  //plots of efficiency with respect to various analysis variables
  if (do_variables) {
	//-------------6.2 table 1 plots (MET and HT dependenece)-------------
	//eff vs MET (real MET)
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_plot_num");
	//eff vs HT (low MET)
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&150<met&&met<=200" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lhtlowmet_plot_den");
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&150<met&&met<=200" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lhtlowmet_plot_num");
	//eff vs HT (high MET)
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met&&met<=300" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lhthighmet_plot_den");
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met&&met<=300" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lhthighmet_plot_num");
	//-------------6.2 table 2 plots (AN variable dependenece- Nj Nb DRmax)-------------
	//eff vs nb
  	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, nb_higgsino, "N_{b}", {-0.5,4.5}),
  	                pass_filters && "njet>=4&&njet<=5&&!low_dphi_met&&nel==1&&200<met" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb_plot_den");
  	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, nb_higgsino, "N_{b}", {-0.5,4.5}),
  	                pass_filters && "njet>=4&&njet<=5&&!low_dphi_met&&nel==1&&200<met" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb_plot_num");
	//eff vs njet
  	pm.Push<Hist1D>(Axis(5, 2.5, 7.5, "njet", "N_{j}", {2.5,7.5}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_njet_plot_den");
  	pm.Push<Hist1D>(Axis(5, 2.5, 7.5, "njet", "N_{j}", {2.5,7.5}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<met" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_njet_plot_num");
	//eff vs drmax
  	pm.Push<Hist1D>(Axis(20, 0., 4., "hig_cand_drmax[0]", "#Delta R_{max}", {0.,4.}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&200<met" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax_plot_den");
  	pm.Push<Hist1D>(Axis(20, 0., 4., "hig_cand_drmax[0]", "#Delta R_{max}", {0.,4.}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&200<met" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax_plot_num");
	//-------------6.2 table 3 plots (AN variable dependenece- <mbb>)-------------
	//eff vs MET (250<= HT< 300, <m> < 100)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300&&hig_cand_am[0]<100" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht250am0_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300&&hig_cand_am[0]<100" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht250am0_plot_num");
	//eff vs MET (250<= HT< 300, 100 <= <m> < 150)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300&&100<=hig_cand_am[0]&&hig_cand_am[0]<150" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht250am100_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300&&100<=hig_cand_am[0]&&hig_cand_am[0]<150" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht250am100_plot_num");
	//eff vs MET (250<= HT< 300, 150 <= <m>)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300&&150<=hig_cand_am[0]" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht250am150_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&250<=ht&&ht<300&&150<=hig_cand_am[0]" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht250am150_plot_num");

	//eff vs MET (300<= HT< 400, <m> < 100)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400&&hig_cand_am[0]<100" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht300am0_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400&&hig_cand_am[0]<100" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht300am0_plot_num");
	//eff vs MET (300<= HT< 400, 100 <= <m> < 150)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400&&100<=hig_cand_am[0]&&hig_cand_am[0]<150" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht300am100_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400&&100<=hig_cand_am[0]&&hig_cand_am[0]<150" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht300am100_plot_num");
	//eff vs MET (300<= HT< 400, 150 <= <m>)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400&&150<=hig_cand_am[0]" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht300am150_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&300<=ht&&ht<400&&150<=hig_cand_am[0]" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht300am150_plot_num");

	//eff vs MET (400<= HT< 600, <m> < 100)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600&&hig_cand_am[0]<100" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht400am0_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600&&hig_cand_am[0]<100" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht400am0_plot_num");
	//eff vs MET (400<= HT< 600, 100 <= <m> < 150)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600&&100<=hig_cand_am[0]&&hig_cand_am[0]<150" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht400am100_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600&&100<=hig_cand_am[0]&&hig_cand_am[0]<150" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht400am100_plot_num");
	//eff vs MET (400<= HT< 600, 150 <= <m>)
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600&&150<=hig_cand_am[0]" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht400am150_plot_den");
  	pm.Push<Hist1D>(Axis(50, 150., 400., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&400<=ht&&ht<600&&150<=hig_cand_am[0]" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metht400am150_plot_num");
	//-------------Misc <mbb> plots not in AN-------------
	//eff vs <m> higgs
  	pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&200<met" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am_plot_den");
  	pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&200<met" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am_plot_num");
	////eff vs <m> higgs 200<HT<300 for HT studies
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&200<ht&&ht<300", procs_met150, all_plot_types).Tag("FixName:trig_raw_amht200_plot_den");
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&200<ht&&ht<300" && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_amht200_plot_num");
	////eff vs <m> higgs 300<HT<400 for HT studies
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&300<ht&&ht<400", procs_met150, all_plot_types).Tag("FixName:trig_raw_amht300_plot_den");
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&300<ht&&ht<400" && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_amht300_plot_num");
	////eff vs <m> higgs 400<HT<500 for HT studies
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&400<ht&&ht<500", procs_met150, all_plot_types).Tag("FixName:trig_raw_amht400_plot_den");
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&400<ht&&ht<500" && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_amht400_plot_num");
	////eff vs <m> higgs 500<HT<600 for HT studies
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&500<ht&&ht<600", procs_met150, all_plot_types).Tag("FixName:trig_raw_amht500_plot_den");
  	//pm.Push<Hist1D>(Axis(20, 0., 300., "hig_cand_am[0]", "#LT m#GT [GeV]", {0.,250.}),
  	//                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=4&&!low_dphi_met&&nel==1&&nmu==0&&150<met&&500<ht&&ht<600" && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am_ht500_plot_num");

	//-------------6.3 table 1 plots (MET and HT dependenece)-------------
	//eff vs MET (QCD MET)
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_plot_num");
	//eff vs HT (low fake MET)
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&150<met&&met<=200" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_htlowmetqcd_plot_den");
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&150<met&&met<=200" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_htlowmetqcd_plot_num");
	//eff vs HT (high fake MET)
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&200<met&&met<=300" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_hthighmetqcd_plot_den");
  	pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "HT [GeV]", {0.,1300.}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&200<met&&met<=300" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_hthighmetqcd_plot_num");
	if (do_controlregions) {
		//-------------6.4 table 1 plots (1l CR plots)-------------
		//eff vs MET (1e region)
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                pass_filters && "njet>=2&&nel==1" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_met1e_plot_den");
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                pass_filters && "njet>=2&&nel==1" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_met1e_plot_num");
		//eff vs el_pt (1el region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Electron p_{T} [GeV]", {20.,170.}),
  		                pass_filters && "njet>=2&&nel==1&&150<met&&met<200" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_elpt1e_plot_den");
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Electron p_{T} [GeV]", {20.,170.}),
  		                pass_filters && "njet>=2&&nel==1&&150<met&&met<200" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_elpt1e_plot_num");
		//eff vs ht (1el region, low MET)
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==1&&150<met&&met<200" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_htlowmet1e_plot_den");
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==1&&150<met&&met<200" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_htlowmet1e_plot_num");
		//eff vs ht (1el region, high MET)
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==1&&200<met&&met<300" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_hthighmet1e_plot_den");
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==1&&200<met&&met<300" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_hthighmet1e_plot_num");
		//eff vs MET (1mu region)
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                pass_filters && "njet>=2&&nmu==1" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_met1mu_plot_den");
  		pm.Push<Hist1D>(Axis(55, 0., 550., "met", "Offline MET [GeV]", {0,1500}),
  		                pass_filters && "njet>=2&&nmu==1" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_met1mu_plot_num");
		//eff vs mu_pt (1mu region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Muon p_{T} [GeV]", {20.,170.}),
  		                pass_filters && "njet>=2&&nmu==1&&150<met&&met<200" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_mupt1mu_plot_den");
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Muon p_{T} [GeV]", {20.,170.}),
  		                pass_filters && "njet>=2&&nmu==1&&150<met&&met<200" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_mupt1mu_plot_num");
		//eff vs ht (1mu region, low MET)
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==1&&150<met&&met<200" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_htlowmet1mu_plot_den");
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==1&&150<met&&met<200" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_htlowmet1mu_plot_num");
		//eff vs ht (1mu region, high MET)
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==1&&200<met&&met<300" && jetht_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_hthighmet1mu_plot_den");
  		pm.Push<Hist1D>(Axis(13, 0., 1300., "ht", "Offline H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==1&&200<met&&met<300" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_hthighmet1mu_plot_num");
		//-------------6.4 table 1 plots (2l CR plots)-------------
		//eff vs el_pt (2el region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_elpt2e_plot_den");
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger) && el_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_elpt2e_plot_num");
		//eff vs ht (2el region)
  		pm.Push<Hist1D>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_ht2e_plot_den");
  		pm.Push<Hist1D>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger) && el_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_ht2e_plot_num");
		//eff vs mu_pt (2mu region)
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_mupt2mu_plot_den");
  		pm.Push<Hist1D>(Axis(30, 20., 170., Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger) && mu_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_mupt2mu_plot_num");
		//eff vs ht (2mu region)
  		pm.Push<Hist1D>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_ht2mu_plot_den");
  		pm.Push<Hist1D>(Axis(10, 0., 1000., "ht", "H_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger) && mu_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_ht2mu_plot_num");
	}
	//-------------Misc MC comparison not in AN-------------
	//MET120 only (nominal)
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0&&nel==1", procs_met150, all_plot_types).Tag("FixName:trig_raw_met120_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0&&nel==1&&(HLT_PFMET120_PFMHT120_IDTight)", procs_met150, all_plot_types).Tag("FixName:trig_raw_met120_plot_num");
	//MET120 MC
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "(pass&&stitch&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0&&nel==1)", procs_met150_mc, all_plot_types).Weight("weight").Tag("FixName:trig_raw_mc_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                "(pass&&stitch&&(HLT_Ele27_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&njet>=3&&!low_dphi_met&&nel==1&&nmu==0&&nel==1&&(HLT_PFMET120_PFMHT120_IDTight))", procs_met150_mc, all_plot_types).Weight("weight").Tag("FixName:trig_raw_mc_plot_num");
	//-------------HT binned plots to determine trigger binning for 0l-------------
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht<=200" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht1_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&ht<=200" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht1_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<=ht&&ht<=350" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht2_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&200<=ht&&ht<=350" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht2_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&350<=ht&&ht<=450" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht3_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&350<=ht&&ht<=450" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht3_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&450<=ht&&ht<=550" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht4_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&450<=ht&&ht<=550" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht4_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&550<=ht&&ht<=650" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht5_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&550<=ht&&ht<=650" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht5_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&650<=ht&&ht<=800" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht6_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&650<=ht&&ht<=800" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht6_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&800<=ht&&ht<=900" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht7_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&800<=ht&&ht<=900" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht7_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&900<=ht&&ht<=1000" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht8_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&900<=ht&&ht<=1000" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht8_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&1000<=ht" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht9_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1&&1000<=ht" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_0lmet_ht9_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&ht<=200" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht1_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&ht<=200" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht1_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&200<=ht&&ht<=350" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht2_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&200<=ht&&ht<=350" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht2_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&350<=ht&&ht<=450" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht3_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&350<=ht&&ht<=450" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht3_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&450<=ht&&ht<=550" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht4_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&450<=ht&&ht<=550" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht4_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&550<=ht&&ht<=650" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht5_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&550<=ht&&ht<=650" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht5_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&650<=ht&&ht<=800" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht6_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&650<=ht&&ht<=800" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht6_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&800<=ht&&ht<=900" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht7_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&800<=ht&&ht<=900" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht7_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&900<=ht&&ht<=1000" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht8_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&900<=ht&&ht<=1000" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht8_plot_num");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&1000<=ht" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht9_plot_den");
  	pm.Push<Hist1D>(Axis(100, 150., 550., "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&1000<=ht" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_metqcd_ht9_plot_num");
  }
  //0l systematics plots
  if (do_systematics) {
	//nj3 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj3_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj3_sys_num");
	//nj4
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet==4&&!low_dphi_met&&nel==1" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj4_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet==4&&!low_dphi_met&&nel==1" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj4_sys_num");
	//nj5
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet==5&&!low_dphi_met&&nel==1" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj5_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet==5&&!low_dphi_met&&nel==1" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nj5_sys_num");
	//nb0
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && nb_higgsino==0., procs_met150, all_plot_types).Tag("FixName:trig_raw_nb0_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && nb_higgsino==0. && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb0_sys_num");
	//nb1
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && nb_higgsino==1., procs_met150, all_plot_types).Tag("FixName:trig_raw_nb1_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && nb_higgsino==1. && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb1_sys_num");
	//nb2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && nb_higgsino>=2., procs_met150, all_plot_types).Tag("FixName:trig_raw_nb2_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && nb_higgsino>=2. && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_nb2_sys_num");
	//drmax < 2.2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]<2.2" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax0_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]<2.2" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax0_sys_num");
	//drmax > 2.2
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]>=2.2" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax22_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&!low_dphi_met&&nel==1&&hig_cand_drmax[0]>=2.2" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_drmax22_sys_num");
	//300<=ht<400 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_ht300to400_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_ht300to400_sys_num");
	//am<140, 300<=ht<400
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]<140" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am0_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]<140" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am0_sys_num");
	//am>140, 300<=ht<400
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]>=140" && el_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am140_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=4&&300<=ht&&ht<400&&!low_dphi_met&&nel==1&&hig_cand_am[0]>=140" && el_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_am140_sys_num");
	//jet_pt>500 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_refel_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "njet>=3&&!low_dphi_met&&nel==1" && el_trigger && high_pt_jet && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_refel_sys_num");
	//jet_pt>500, ref trigger jet500
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_refjet_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_refjet_sys_num");
	//jet_pt>500 nel=1, ref trigger jet500 (nominal)
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_1e_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nel==1" && high_pt_jet && metnomu_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_1e_sys_num");
	//jet_pt>500 nmu=1, ref trigger jet500
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nmu==1" && high_pt_jet, procs_met150, all_plot_types).Tag("FixName:trig_raw_1mu_sys_den");
  	pm.Push<Hist1D>(Axis(sys_met_bins, "met", "MET [GeV]", {150,1500}),
  	                pass_filters && "(met/mht<2)&&(met/met_calo<2)&&HLT_PFJet500&&njet>=3&&!low_dphi_met&&nmu==1" && high_pt_jet && metnomu_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_1mu_sys_num");
  }
  // trigger efficiency plots and application file
  if (do_efficiency) {
  	pm.Push<Hist2D>(Axis(true_met_bins, "met", "MET", {}),
  	                Axis(ht_bins, "ht", "HT", {}),
  	                pass_filters && "!low_dphi_met&&njet>=4&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && el_trigger, procs_met150, twodim_plotopts).Tag("FixName:trig_raw_0l_eff_den");
  	pm.Push<Hist2D>(Axis(true_met_bins, "met", "MET", {}),
  	                Axis(ht_bins, "ht", "HT", {}),
  	                pass_filters && "!low_dphi_met&&njet>=4&&nel==1" && (Higfuncs::lead_signal_lepton_pt<35) && el_trigger && met_trigger, procs_met150, twodim_plotopts).Tag("FixName:trig_raw_0l_eff_num");
	//QCD HT bins
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht0to350, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&ht<350" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht0to350_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht0to350, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&ht<350" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht0to350_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht350to450, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&350<=ht&&ht<450" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht350to450_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht350to450, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&350<=ht&&ht<450" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht350to450_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht450to550, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&450<=ht&&ht<550" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht450to550_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht450to550, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&450<=ht&&ht<550" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht450to550_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht550to650, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&550<=ht&&ht<650" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht550to650_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins_ht550to650, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&550<=ht&&ht<650" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht550to650_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&650<=ht&&ht<800" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht650to800_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&650<=ht&&ht<800" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht650to800_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&800<=ht&&ht<1000" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht800to1000_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&800<=ht&&ht<1000" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht800to1000_num");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&1000<=ht" && jetht_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht1000toInf_den");
  	pm.Push<Hist1D>(Axis(fake_met_bins, "met", "MET [GeV]", {}),
  	                pass_filters && "low_dphi_met&&nvlep==0&&1000<=ht" && jetht_trigger && met_trigger, procs_met150, all_plot_types).Tag("FixName:trig_raw_qcd_eff_ht1000toInf_num");
  	//pm.Push<Hist2D>(Axis(fake_met_bins, "met", "MET", {}),
  	//                Axis(ht_fake_met_bins, "ht", "HT", {}),
  	//                pass_filters && "low_dphi_met&&nvlep==0" && ht_trigger, procs_met150, twodim_plotopts).Tag("FixName:trig_raw_qcd_eff_den");
  	//pm.Push<Hist2D>(Axis(fake_met_bins, "met", "MET", {}),
  	//                Axis(ht_fake_met_bins, "ht", "HT", {}),
  	//                pass_filters && "low_dphi_met&&nvlep==0" && ht_trigger && met_trigger, procs_met150, twodim_plotopts).Tag("FixName:trig_raw_qcd_eff_num");
	if (do_controlregions) {
		//1e
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                pass_filters && "njet>=2&&nel==1&&0<=ht&&ht<400" && jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht0to400_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                pass_filters && "njet>=2&&nel==1&&0<=ht&&ht<400" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht0to400_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                pass_filters && "njet>=2&&nel==1&&400<=ht&&ht<600" && jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht400to600_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                pass_filters && "njet>=2&&nel==1&&400<=ht&&ht<600" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht400to600_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                pass_filters && "njet>=2&&nel==1&&600<=ht" && jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht600toInf_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(el_pt_bins, Higfuncs::lead_signal_electron_pt, "Electron pt", {}),
  		                pass_filters && "njet>=2&&nel==1&&600<=ht" && jetht_trigger && (met_trigger||el_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1el_ht600toInf_eff_num");
		//1mu
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                pass_filters && "njet>=2&&nmu==1&&0<=ht&&ht<=400" && jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht0to400_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                pass_filters && "njet>=2&&nmu==1&&0<=ht&&ht<=400" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht0to400_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                pass_filters && "njet>=2&&nmu==1&&400<=ht&&ht<=600" && jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht400to600_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                pass_filters && "njet>=2&&nmu==1&&400<=ht&&ht<=600" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht400to600_eff_num");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                pass_filters && "njet>=2&&nmu==1&&600<=ht" && jetht_trigger, procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht600toInf_eff_den");
  		pm.Push<Hist2D>(Axis(twodim_met_bins, "met", "MET", {}),
  		                Axis(mu_pt_bins, Higfuncs::lead_signal_muon_pt, "Muon pt", {}),
  		                pass_filters && "njet>=2&&nmu==1&&600<=ht" && jetht_trigger && (met_trigger||mu_trigger), procs_1l2j, twodim_plotopts).Tag("FixName:trig_raw_1mu_ht600toInf_eff_num");
		//2l
  		pm.Push<Hist1D>(Axis(twoel_pt_bins, Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2el_eff_den");
  		pm.Push<Hist1D>(Axis(twoel_pt_bins, Higfuncs::lead_signal_electron_pt, "Offline Max Electron p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nel==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger) && el_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2el_eff_num");
  		pm.Push<Hist1D>(Axis(twomu_pt_bins, Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger), procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2mu_eff_den");
  		pm.Push<Hist1D>(Axis(twomu_pt_bins, Higfuncs::lead_signal_muon_pt, "Offline Max Muon p_{T} [GeV]", {}),
  		                pass_filters && "njet>=2&&nmu==2&&(80<ll_m[0]&&ll_m[0]<100)" && (met_trigger||jetht_trigger) && mu_trigger, procs_1l2j, all_plot_types).Tag("FixName:trig_raw_2mu_eff_num");
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
	//-------------6.2 table 1 plots (MET and HT dependenece)-------------
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_"+year_string+".pdf");
	pm_idx += 2;
	
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 150<MET#leq 200, 35.9 fb^{-1} (13 TeV); Offline HT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_htlowmet");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 150 GeV #leq p_{T}^{miss}<200 GeV","Offline H_{T}","Efficiency "+str_met_trig,"GeV","htlowmet_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline HT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_hthighmet");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 200 GeV #leq p_{T}^{miss}<300 GeV","Offline H_{T}","Efficiency "+str_met_trig,"GeV","hthighmet_"+year_string+".pdf");
	pm_idx += 2;

	//-------------6.2 table 2 plots (AN variable dependenece- Nj Nb DRmax)-------------
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline N_{b}; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nb");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV","Offline N_{b}","Efficiency "+str_met_trig,"","nb_nj4_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline N_{j}; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_nj");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV","Offline N_{j}","Efficiency "+str_met_trig,"","nj_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline #Delta R_{max}; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcanddrmax");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV","Offline #Delta R_{max}","Efficiency "+str_met_trig,"","higcanddrmax_"+year_string+".pdf");
	pm_idx += 2;

	//-------------6.2 table 3 plots (AN variable dependenece- <mbb>)-------------
	make_ambb_plots(&pm, pm_idx, pm_idx+1, pm_idx+2, pm_idx+3, pm_idx+4, pm_idx+5, pro_met150, str_lumi, "250 #leq H_{T} < 300 GeV", "ht250to300_"+year_string);
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 250#leq H_{T} < 300 GeV #LT m_{bb}#GT < 100 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb0to100");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 250 #leq H_{T} < 300 GeV, #LT m_{bb} #GT < 100 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb0to100ht250to300_"+year_string+".pdf");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 250#leq H_{T} < 300 GeV 100 #leq #LT m_{bb}#GT < 150 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb100to150");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 250 #leq H_{T} < 300 GeV, 100 #leq #LT m_{bb} #GT < 150 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb100to150ht250to300_"+year_string+".pdf");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 250#leq H_{T} < 300 GeV 150 GeV #leq #LT m_{bb}#GT, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb150toInf");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 250 #leq H_{T} < 300 GeV, 150 GeV #leq #LT m_{bb} #GT","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb150toinfht250to300_"+year_string+".pdf");
	pm_idx += 2;
	
	make_ambb_plots(&pm, pm_idx, pm_idx+1, pm_idx+2, pm_idx+3, pm_idx+4, pm_idx+5, pro_met150, str_lumi, "300 #leq H_{T} < 400 GeV", "ht300to400_"+year_string);
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 300#leq H_{T} < 400 GeV #LT m_{bb}#GT < 100 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb0to100");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 300 #leq H_{T} < 400 GeV, #LT m_{bb} #GT < 100 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb0to100ht300to400_"+year_string+".pdf");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 300#leq H_{T} < 400 GeV 100 #leq #LT m_{bb}#GT < 150 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb100to150");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 300 #leq H_{T} < 400 GeV, 100 #leq #LT m_{bb} #GT < 150 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb100to150ht300to400_"+year_string+".pdf");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 300#leq H_{T} < 400 GeV 150 GeV #leq #LT m_{bb}#GT, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb150toInf");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 300 #leq H_{T} < 400 GeV, 150 GeV #leq #LT m_{bb} #GT","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb150toinfht300to400_"+year_string+".pdf");
	pm_idx += 2;

	make_ambb_plots(&pm, pm_idx, pm_idx+1, pm_idx+2, pm_idx+3, pm_idx+4, pm_idx+5, pro_met150, str_lumi, "400 #leq H_{T} < 600 GeV", "ht400to600_"+year_string);
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 400 #leq H_{T} < 600 GeV #LT m_{bb}#GT < 100 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb0to100");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 400 #leq H_{T} < 600 GeV, #LT m_{bb} #GT < 100 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb0to100ht400to600_"+year_string+".pdf");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 400 #leq H_{T} < 600 GeV 100 #leq #LT m_{bb}#GT < 150 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb100to150");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 400 #leq H_{T} < 600 GeV, 100 #leq #LT m_{bb} #GT < 150 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb100to150ht400to600_"+year_string+".pdf");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 4 high #Delta#phi N_{e}=1 400 #leq H_{T} < 600 GeV 150 GeV #leq #LT m_{bb}#GT, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_metambb150toInf");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, 400 #leq H_{T} < 600 GeV, 150 GeV #leq #LT m_{bb} #GT","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","met_ambb150toinfht400to600_"+year_string+".pdf");
	pm_idx += 2;

	//-------------Misc <mbb> plots not in AN-------------
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 200<MET, 35.9 fb^{-1} (13 TeV); Offline #LT m#RT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcandam");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV","Offline #LT m_{bb} #GT","Efficiency "+str_met_trig,"GeV","higcandam_"+year_string+".pdf");
	pm_idx += 2;

  	//generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 150<MET 200<HT<300, 35.9 fb^{-1} (13 TeV); Offline #LT m#RT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcandamht200300");
	//make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV, 200 #leq H_{T} < 300 GeV","Offline #LT m_{bb} #GT","Efficiency "+str_met_trig,"GeV","higcandam_ht200to300_"+year_string+".pdf");
	//pm_idx += 2;
	//
  	//generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 150<MET 300<HT<400, 35.9 fb^{-1} (13 TeV); Offline #LT m#RT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcandamht300400");
	//make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV, 300 #leq H_{T} < 400 GeV","Offline #LT m_{bb} #GT","Efficiency "+str_met_trig,"GeV","higcandam_ht300to400_"+year_string+".pdf");
	//pm_idx += 2;

  	//generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 150<MET 400<HT<500, 35.9 fb^{-1} (13 TeV); Offline #LT m#RT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcandamht400500");
	//make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV, 400 #leq H_{T} < 500 GeV","Offline #LT m_{bb} #GT","Efficiency "+str_met_trig,"GeV","higcandam_ht400to500_"+year_string+".pdf");
	//pm_idx += 2;

  	//generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1 150<MET 500<HT<600, 35.9 fb^{-1} (13 TeV); Offline #LT m#RT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_higcandamht500600");
	//make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 4, high #Delta #phi, N_{e}=1, p_{T}^{miss} #geq 200 GeV, 500 #leq H_{T} < 600 GeV","Offline #LT m_{bb} #GT","Efficiency "+str_met_trig,"GeV","higcandam_ht500to600_"+year_string+".pdf");
	//pm_idx += 2;

	//-------------6.3 table 1 plots (QCD CR MET and HT dependenece)-------------
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=1 150<MET#leq 200, 35.9 fb^{-1} (13 TeV); Offline HT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_htfakelowmet");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 150 #leq p_{T}^{miss} < 200 GeV","Offline H_{T}","Efficiency "+str_met_trig,"GeV","htfakelowmet_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0 200<MET, 35.9 fb^{-1} (13 TeV); Offline HT [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_htfakehighmet");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 200 #leq p_{T}^{miss} < 300 GeV","Offline H_{T}","Efficiency "+str_met_trig,"GeV","htfakehighmet_"+year_string+".pdf");
	pm_idx += 2;

	if (do_controlregions) {
		//-------------6.4 table 1 plots (1l CR plots)-------------
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||Ele(27_WPTight||35_WPTight||115)]","hist_elmet");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{e}=1","Offline p_{T}^{miss}","Efficiency "+str_met_trig+"|"+str_el_trig,"GeV","elmet_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline Electron Pt [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||Ele(27_WPTight||35_WPTight||115)]","hist_elpt");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{e}=1, 150 #leq p_{T}^{miss} < 200 GeV","Offline Electron p_{T}","Efficiency "+str_met_trig+"|"+str_el_trig,"GeV","elpt_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline H_{T} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||Ele(27_WPTight||35_WPTight||115)]","hist_elhtlowmet");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{e}=1, 150 #leq p_{T}^{miss} < 200 GeV","Offline H_{T}","Efficiency "+str_met_trig+"|"+str_el_trig,"GeV","el_htlowmet_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline H_{T} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||Ele(27_WPTight||35_WPTight||115)]","hist_elhthighmet");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{e}=1, 200 #leq p_{T}^{miss} < 300 GeV","Offline H_{T}","Efficiency "+str_met_trig+"|"+str_el_trig,"GeV","el_hthighmet_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||(IsoMu(24||27)||Mu50)]","hist_mumet");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{#mu}=1","Offline p_{T}^{miss}","Efficiency "+str_met_trig+"|"+str_mu_trig,"GeV","mumet_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline Muon Pt [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||(IsoMu(24||27)||Mu50)]","hist_mupt");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{#mu}=1, 150 #leq p_{T}^{miss} < 200 GeV","Offline Muon p_{T}","Efficiency "+str_met_trig+"|"+str_mu_trig,"GeV","mupt_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline H_{T} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||(IsoMu(24||27)||Mu50)]","hist_muhtlowmet");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{#mu}=1, 150 #leq p_{T}^{miss} < 200 GeV","Offline H_{T}","Efficiency "+str_met_trig+"|"+str_mu_trig,"GeV","mu_htlowmet_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 2 N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline H_{T} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)||(IsoMu(24||27)||Mu50)]","hist_muhtlowmet");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{JetHT}, N_{j}#geq 2, N_{#mu}=1, 200 #leq p_{T}^{miss} < 300 GeV","Offline H_{T}","Efficiency "+str_met_trig+"|"+str_mu_trig,"GeV","mu_hthighmet_"+year_string+".pdf");
		pm_idx += 2;

		//-------------6.4 table 1 plots (2l CR plots)-------------
  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{e}=2 80<m_{ee}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Electron Pt [GeV]; Efficiency [Ele(27_WPTight||35_WPTight||115)]","hist_elel_show");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{e}=2, 80 < m_{ee} #leq 100 GeV","Max Offline Electron p_{T}","Efficiency "+str_el_trig,"GeV","elel_show_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{e}=2 80<m_{ee}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline H_{T} [GeV]; Efficiency [Ele(27_WPTight||35_WPTight||115)]","hist_elel_ht");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{e}=2, 80 < m_{ee} #leq 100 GeV","Offline H_{T}","Efficiency "+str_el_trig,"GeV","elel_ht_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{#mu}=2 80<m_{#mu#mu}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline Max Muon Pt [GeV]; Efficiency [IsoMu(24||27)||Mu50]","hist_mumu_show");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{#mu}=2, 80 < m_{#mu#mu} #leq 100 GeV","Max Offline Muon p_{T}","Efficiency "+str_mu_trig,"GeV","mumu_show_"+year_string+".pdf");
		pm_idx += 2;

  		generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_1l2j, "Trigger Efficiency, baseline: HLT_PFJet500||HLT_MET[NoMu](110||120||120_HT60) N_{j}#geq 2 N_{#mu}=2 80<m_{#mu#mu}<100 GeV, 35.9 fb^{-1} (13 TeV); Offline H_{T} [GeV]; Efficiency [IsoMu(24||27)||Mu50]","hist_mumu_ht");
		make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_1l2j,str_lumi,"Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{#mu}=2, 80 < m_{#mu#mu} #leq 100 GeV","H_{T}","Efficiency "+str_mu_trig,"GeV","mumu_ht_"+year_string+".pdf");
		pm_idx += 2;

	}
	//-------------Misc MC comparison not in AN-------------
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET120]","hist_datamet120");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1","Offline p_{T}^{miss}","Efficiency MET[NoMu]120","GeV","datamet120_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150_mc, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, Summer16 t#bar{t} 1l; Offline E_{T}^{miss} [GeV]; Efficiency [MET120]","hist_mcmet120");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150_mc,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1","Offline p_{T}^{miss}","Efficiency MET[NoMu]120","GeV","mcmet120_"+year_string+".pdf",false,true);
	pm_idx += 2;

	//-------------HT binned plots to determine trigger binning for 0l-------------
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, H_{T}<200 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht0to200");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, H_{T}<200 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht0to200_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 200<H_{T}<350 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht200to350");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 200 #leq H_{T}<350 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht200to350_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 350<H_{T}<450 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht350to450");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 350 #leq H_{T}<450 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht350to450_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 450<H_{T}<550 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht450to550");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 450 #leq H_{T}<550 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht450to550_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 550<H_{T}<650 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht550to650");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 550 #leq H_{T}<650 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht550to650_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 650<H_{T}<800 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht650to800");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 650 #leq H_{T}<800 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht650to800_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 800<H_{T}<900 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht800to900");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 800 #leq H_{T}<900 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht800to900_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 900<H_{T}<1000 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht900to1000");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, 900 #leq H_{T}<1000 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht900to1000_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, H_{T}>1000 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_realmet_ht1000toInf");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_realmet_ref_trig+", N_{j}#geq 3, high #Delta #phi, N_{e}=1, H_{T} #geq 1000 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","realmet_ht1000toInf_"+year_string+".pdf");
	pm_idx += 2;

	//QCD
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, H_{T}<200 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht0to200");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, H_{T}<200 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht0to200_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 200<H_{T}<350 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht200to350");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 200 #leq H_{T}<350 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht200to350_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 350<H_{T}<450 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht350to450");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 350 #leq H_{T}<450 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht350to450_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 450<H_{T}<550 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht450to550");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 450 #leq H_{T}<550 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht450to550_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 550<H_{T}<650 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht550to650");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 550 #leq H_{T}<650 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht550to650_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 650<H_{T}<800 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht650to800");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 650 #leq H_{T}<800 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht650to800_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 800<H_{T}<900 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht800to900");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 800 #leq H_{T}<900 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht800to900_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, 900<H_{T}<1000 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht900to1000");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, 900 #leq H_{T}<1000 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht900to1000_"+year_string+".pdf");
	pm_idx += 2;

  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 low #Delta#phi N_{vl}=0, H_{T}>1000 GeV, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_fakemet_ht1000toInf");
	make_efficiency_plot_wrapper(&pm,pm_idx,pm_idx+1,pro_met150,str_lumi,"Denom: "+str_fakemet_ref_trig+", low #Delta #phi, N_{l veto}=0, H_{T} #geq 1000 GeV","Offline p_{T}^{miss}","Efficiency "+str_met_trig,"GeV","fakemet_ht1000toInf_"+year_string+".pdf");
	pm_idx += 2;
  }
  std::vector<double> systematics;
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
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 3 p_{Tjet}>500 GeV, high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_1e");
	pm_idx += 2;
  	generate_ratio_plots(&pm, pm_idx, pm_idx+1, pro_met150, "Trigger Efficiency, baseline: HLT_PFJet500 N_{j}#geq 3 p_{Tjet}>500 GeV, high #Delta#phi N_{#mu}=1, 35.9 fb^{-1} (13 TeV); Offline E_{T}^{miss} [GeV]; Efficiency [MET[NoMu](110||120||120_HT60)]","hist_1mu");
	pm_idx += 2;
  }
  if (do_efficiency) {
	trig_postprocess_efficiencies(&pm, pm_idx, pro_met150, pro_1l2j, year_string, "pdf", systematics);
  	generate_2d_efficiencies(&pm, pm_idx, pm_idx+1, "[MET[NoMu](110|120||120_HT60)] Trigger Efficiency, baseline: HLT_Ele27 N_{j}#geq 3 high #Delta#phi N_{e}=1, 35.9 fb^{-1} (13 TeV); MET [GeV], HT [GeV]", "hist_realmetht");
	pm_idx += 2;

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
       //nel1
       Hist1D * den_h_nel1 = static_cast<Hist1D*>(pm->Figures()[first_index+26].get());
       Hist1D::SingleHist1D* den_sh_nel1 = static_cast<Hist1D::SingleHist1D*>(den_h_nel1->GetComponent(proc.get()));
       TH1D den_nel1 = den_sh_nel1->scaled_hist_;
       Hist1D * num_h_nel1 = static_cast<Hist1D*>(pm->Figures()[first_index+27].get());
       Hist1D::SingleHist1D* num_sh_nel1 = static_cast<Hist1D::SingleHist1D*>(num_h_nel1->GetComponent(proc.get()));
       TH1D num_nel1 = num_sh_nel1->scaled_hist_;
       TGraphAsymmErrors* plot_ratio_nel1 = new TGraphAsymmErrors(&num_nel1,&den_nel1,"cp");
       //nel0
       Hist1D * den_h_nel0 = static_cast<Hist1D*>(pm->Figures()[first_index+28].get());
       Hist1D::SingleHist1D* den_sh_nel0 = static_cast<Hist1D::SingleHist1D*>(den_h_nel0->GetComponent(proc.get()));
       TH1D den_nel0 = den_sh_nel0->scaled_hist_;
       Hist1D * num_h_nel0 = static_cast<Hist1D*>(pm->Figures()[first_index+29].get());
       Hist1D::SingleHist1D* num_sh_nel0 = static_cast<Hist1D::SingleHist1D*>(num_h_nel0->GetComponent(proc.get()));
       TH1D num_nel0 = num_sh_nel0->scaled_hist_;
       TGraphAsymmErrors* plot_ratio_nel0 = new TGraphAsymmErrors(&num_nel0,&den_nel0,"cp");
	
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
	table_file.open("tables/"+out_filename+"_trig_sys.tex");
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
	table_file << "\\multicolumn{8}{c}{$N_\\text{l}$ variations (Jet500 trigger, Max $p_\\text{Tjet}>500$~GeV)}\\\\ \\hline \n";
	table_file << make_systematic_table_row(plot_ratio_nel1,"$N_\\text{e}=1$ (Nominal)") << "\n";
	nel_table.push_back(get_systematic_table_values(plot_ratio_nel1));
	table_file << make_systematic_table_row(plot_ratio_nel0,"$N_\\text{\\mu}=1$") << " \\hline \n";
	nel_table.push_back(get_systematic_table_values(plot_ratio_nel0));
	nel_row = get_systematic_table_variations(nel_table);
	greater_variations_rows.push_back(nel_row);
	table_file << make_systematic_table_row(nel_row) << "\n";
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
	//	table_file << "\\multicolumn{8}{c}{Data compared to MC (MET120\\_HT60 trigger only)}\\\\ \\hline \n";
	//	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio")),"Data (nominal)") << "\n";
	//	datamc_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio"))));
	//	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio")),"MC") << "\\hline \n";
	//	datamc_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio"))));
	//	datamc_row = get_systematic_table_variations(datamc_table);
	//	systematics_rows.push_back(datamc_row);
	//	table_file << make_systematic_table_row(datamc_row) << "\n";
	//}
	std::vector<double> ret = add_variations_quadrature(systematics_rows);
	table_file << make_systematic_table_row(ret, true) << "\n";
	table_file << "\\end{tabular} \n";
        table_file.close();
        std::cout << "Wrote tables/higgsino/"+out_filename+"_"+std::to_string(year)+"_trig_sys.tex" << std::endl;
	return ret;
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
void trig_postprocess_efficiencies(PlotMaker* pm, int first_index, shared_ptr<Process> met_proc, shared_ptr<Process> lep_proc, std::string year, std::string img_file_extension, std::vector<double> systematics) {
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
  draw_2d_efficiency_graphs(pm, first_index, first_index+1, ht_bins, true_met_bins, "ht", "H_{T}", "b.ht()", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{SingleElectron}, N_{j}#geq 3, high #Delta #phi, N_{e}=1", lumi_string, met_trigger_string, "truemet_effhtbin", apply_effs_file, "get_0l_trigeff", year, img_file_extension, systematics);
  //fake met histograms
  draw_2d_efficiency_graphs_variable_binning(pm, {first_index+2,first_index+4,first_index+6,first_index+8,first_index+10,first_index+12,first_index+14}, {first_index+3,first_index+5,first_index+7,first_index+9,first_index+11,first_index+13,first_index+15}, met_proc, ht_fake_met_bins, "ht", "H_{T}", "b.ht()", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{JetHT}, low #Delta #phi, N_{l veto}=0", lumi_string, met_trigger_string, "fakemet_effhtbin", apply_effs_file, "get_0l_fakemet_trigeff", year, img_file_extension);
  if (do_controlregions) {
    //1/2l CR graphs
    draw_3d_efficiency_graphs(pm, {first_index+16, first_index+18, first_index+20}, {first_index+17, first_index+19, first_index+21}, singlelep_ht_bins, el_pt_bins, twodim_met_bins, "ht", "H_{T}", "b.ht()", "el_pt", "p_{Te}", "HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig())", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{JetHT}, N_{j}#geq 2, N_{e}=1", lumi_string, met_trigger_string+"|"+el_trigger_string, "elmet_effelptbin", apply_effs_file, "get_1el_trigeff", year, img_file_extension);
    draw_3d_efficiency_graphs(pm, {first_index+22, first_index+24, first_index+26}, {first_index+23, first_index+25, first_index+27}, singlelep_ht_bins, mu_pt_bins, twodim_met_bins, "ht", "H_{T}", "b.ht()", "mu_pt", "p_{T#mu}", "HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig())", "met", "p_{T}^{miss}", "b.met()", "Denom: #it{JetHT}, N_{j}#geq 2, N_{#mu}=1", lumi_string, met_trigger_string+"|"+mu_trigger_string, "mumet_effmuptbin", apply_effs_file, "get_1mu_trigeff", year, img_file_extension);
    draw_1d_efficiency_graphs(pm, first_index+28, first_index+29, lep_proc, twoel_pt_bins, "el_pt", "Max p_{Te}", "HigUtilities::signal_lepton_pt(b.el_pt(),b.el_sig())","Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{e}=2",lumi_string,el_trigger_string,"elel_eff",apply_effs_file,"get_2el_trigeff",year,img_file_extension);
    draw_1d_efficiency_graphs(pm, first_index+30, first_index+31, lep_proc, twomu_pt_bins, "mu_pt", "Max p_{T#mu}", "HigUtilities::signal_lepton_pt(b.mu_pt(),b.mu_sig())","Denom: #it{MET|JetHT}, N_{j}#geq 2, N_{#mu}=2",lumi_string,mu_trigger_string,"mumu_eff",apply_effs_file,"get_2mu_trigeff",year,img_file_extension);
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

void draw_2d_efficiency_graphs_variable_binning(PlotMaker* pm, std::vector<int> den_idx, std::vector<int> num_idx, shared_ptr<Process> proc, std::vector<double> y_bins, std::string y_var_name, std::string y_title_name, std::string y_expression, std::string x_var_name, std::string x_title_name, std::string x_expression, std::string cuts_string, std::string lumi_string, std::string new_trigger_string, std::string plot_name, ofstream &out_file, std::string func_name, std::string year, std::string img_file_extension) {
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
      {"noplots", no_argument, 0, 0},
      {"nosystematics", no_argument, 0, 0},
      {"noefficiency", no_argument, 0, 0},
      {"nocr", no_argument, 0, 0},
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
      } else if (optname == "noplots") {
        do_variables = false;
      } else if (optname == "nosystematics") {
        do_systematics = false;
      } else if (optname == "noefficiency") {
        do_efficiency = false;
      } else if (optname == "nocr") {
        do_controlregions = false;
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
